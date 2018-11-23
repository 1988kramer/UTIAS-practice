addpath ../common;

deltaT = .02;
alphas = [.11 .01 .18 .08 0 0]; % motion model noise parameters

% measurement model noise parameters
Q_t = [11.7   0;
       0    0.18];
   
n_robots = 1;
robot_num = 1;
[Barcodes, Landmark_Groundtruth, Robots] = loadMRCLAMdataSet(n_robots);
[Robots, timesteps] = sampleMRCLAMdataSet(Robots, deltaT);

% add pose estimate matrix to Robots
Robots{robot_num}.Est = zeros(size(Robots{robot_num}.G,1), 4);

start = 600;
t = Robots{robot_num}.G(start, 1);

% initialize state mean
stateMean = [Robots{robot_num}.G(start,2);
            Robots{robot_num}.G(start,3);
            Robots{robot_num}.G(start,4)];
        
stateCov = zeros(3, 3);
stateCov(1:3,1:3) = 0.001;

measurementIndex = 1;
meas_count = 0;
n_landmarks = 0;

% set up map between barcodes and landmark IDs
codeDict = containers.Map(Barcodes(:,2),Barcodes(:,1));

% loop through all odometry and measurement samples
% updating the robot's pose estimate with each step
% reference table 10.1 in Probabilistic Robotics
for i = start:size(Robots{robot_num}.G, 1)
    
    % update time
    t = Robots{robot_num}.G(i, 1);
    
    % update movement vector
    u_t = [Robots{robot_num}.O(i, 2); Robots{robot_num}.O(i, 3)];
    % update robot bearing
    theta = stateMean(3, 1);

    rot = deltaT * u_t(2);
    halfRot = rot / 2;
    trans = u_t(1) * deltaT;
    
    % calculate pose update from odometry
    poseUpdate = [trans * cos(theta + halfRot);
                  trans * sin(theta + halfRot);
                  rot];
              
    % calculate updated state mean
    F_x = [eye(3) zeros(3, size(stateMean, 1) - 3)];
    stateMeanBar = stateMean + F_x' * poseUpdate;
    stateMeanBar(3) = conBear(stateMeanBar(3));
    
    % calculate movement jacobian
    g_t = [0 0 trans * -sin(theta + halfRot);
           0 0 trans * cos(theta + halfRot);
           0 0 0];
    G_t = eye(size(stateMean, 1)) + F_x' * g_t * F_x;
    
    % calculate motion covariance in control space
    M_t = [(alphas(1) * abs(u_t(1)) + alphas(2) * abs(u_t(2)))^2 0;
           0 (alphas(3) * abs(u_t(1)) + alphas(4) * abs(u_t(2)))^2];
       
    % calculate Jacobian to transform motion covariance to state space
    V_t = [cos(theta + halfRot) -0.5 * sin(theta + halfRot);
           sin(theta + halfRot) 0.5 * cos(theta + halfRot);
           0 1];
    
    % update state covariance
    R_t = V_t * M_t * V_t';
    stateCovBar = (G_t * stateCov * G_t') + (F_x' * R_t * F_x);
    
    % get measurements
    [z, measurementIndex] = getObservations(Robots, robot_num, t, measurementIndex, codeDict);
    
    % if features are observed
    % loop over all features and compute Kalman gain
    if z(3,1) > 0
        for k = 1:size(z, 2) % loop over every measurement

            predZ   = zeros(2, n_landmarks+1);
            predPsi = zeros(n_landmarks+1, 2, 2);
            predH   = zeros(n_landmarks+1, 2, 2*(n_landmarks+1)+3);
            pi_k    = zeros(n_landmarks+1, 1);
            
            % create temporary new landmark at observed position
            temp_mark = [stateMeanBar(1) + z(1,k) * cos(z(2,k) + stateMeanBar(3));
                         stateMeanBar(2) + z(1,k) * sin(z(2,k) + stateMeanBar(3))]; 
            stateMeanTemp = [stateMeanBar;
                             temp_mark];
            stateCovTemp = [stateCovBar zeros(size(stateCovBar,1),2);
                            zeros(2,size(stateCovBar,2) + 2)];
            % initialize covariance for new landmark proportional
            % to range measurement squared
            for ii = (size(stateCovTemp,1)-1):(size(stateCovTemp,1))
                stateCovTemp(ii,ii) = (z(1,k)^2)/130;
            end
            
            % loop over all landmarks  (including the temp landmark) and 
            % compute likelihood of correspondence with the new landmark
            % NOTE: could improve by caching predicted observations when
            %       more than 1 observation occurs at the same timestep
            max_j = 0;
            min_pi = 10*ones(2,1);
            for j = 1:n_landmarks+1

                delta = [stateMeanTemp(2*j+2) - stateMeanTemp(1);
                         stateMeanTemp(2*j+3) - stateMeanTemp(2)];
                       
                q = delta' * delta;
                r = sqrt(q);
                
                predZ(:,j) = [r;
                              conBear(atan2(delta(2), delta(1)) - stateMeanTemp(3))];
                F_xj = [eye(3)     zeros(3,2*j-2) zeros(3,2) zeros(3,2*(n_landmarks+1) - 2*j);
                        zeros(2,3)   zeros(2,2*j-2) eye(2)   zeros(2,2*(n_landmarks+1) - 2*j)];
                
                
                h_t = [-delta(1)/r -delta(2)/r  0   delta(1)/r delta(2)/r;
                       delta(2)/q    -delta(1)/q    -1  -delta(2)/q  delta(1)/q];
                predH(j,:,:) = h_t * F_xj;
                predPsi(j,:,:) = squeeze(predH(j,:,:)) * stateCovTemp * ...
                                 squeeze(predH(j,:,:))' + Q_t;
                if j <= n_landmarks
                    pi_k(j) = ((z(1:2,k)-predZ(:,j))'...
                              /(squeeze(predPsi(j,:,:))))...
                              *(z(1:2,k)-predZ(:,j));
                else
                    pi_k(j) = 0.84; % alpha: min mahalanobis distance to
                                    %        add landmark to map
                end
                % track best two associations
                if pi_k(j) < min_pi(1)
                    min_pi(2) = min_pi(1);
                    max_j = j;
                    min_pi(1) = pi_k(j);
                end
            end
            
            H = squeeze(predH(max_j,:,:));
            
            % best association must be significantly better than
            % the second best, otherwise the measurement is 
            % thrown out
            if (min_pi(2) / min_pi(1) > 1.6)
                meas_count = meas_count + 1;
                % if a landmark is added to the map expand the 
                % state and covariance matrices
                if max_j > n_landmarks
                    stateMeanBar = stateMeanTemp;
                    stateCovBar = stateCovTemp;
                    n_landmarks = n_landmarks + 1;
                % if measurement is associated to an existing landmark
                % truncate h matrix to prevent dim mismatch
                else
                    H = H(:,1:n_landmarks*2 + 3);

                    K = stateCovBar * H' / (squeeze(predPsi(max_j,:,:))); 

                    stateMeanBar = stateMeanBar + K * ...
                                    (z(1:2,k) - predZ(:,max_j));
                    stateMeanBar(3) = conBear(stateMeanBar(3));
                    stateCovBar = (eye(size(stateCovBar)) - (K * H)) * stateCovBar;
                end
            end
        end
    end
    
    % update state mean and covariance
    stateMean = stateMeanBar;
    stateCov = stateCovBar;
    
    % add new pose mean to estimated poses
    Robots{robot_num}.Est(i,:) = [t stateMean(1) stateMean(2) stateMean(3)];
end

p_loss = path_loss(Robots, robot_num, start);
disp(p_loss);
Robots{robot_num}.Map = stateMean;

%animateMRCLAMdataSet(Robots, Landmark_Groundtruth, timesteps, deltaT);