addpath ../common;

deltaT = .02;
alphas = [.1 .01 .18 .08 0 0]; % motion model noise parameters

% measurement model noise parameters
Q_t = [11.8 0    0;
       0    .18  0;
       0    0    1];
   
n_robots = 1;
robot_num = 1;
n_landmarks = 15;
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

% set up map between barcodes and landmark IDs
codeDict = containers.Map(Barcodes(:,2),Barcodes(:,1));

% loop through all odometry and measurement samples
% updating the robot's pose estimate with each step
% reference table 10.1 in Probabilistic Robotics
for i = start:size(Robots{robot_num}.G, 1)
    
    % update time
    t = Robots{robot_num}.G(i, 1);
    
    % update number of landmarks
    n_landmarks = size(stateMean, 1) - 3;
    
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
    if z(3,1) > 1
        for k = 1:size(z, 2) % loop over every observed landmark

            predZ   = zeros(n_landmarks+1, 1, 3);
            predPsi = zeros(n_landmarks+1, 3, 3);
            predH   = zeros(n_landmarks+1, 3, 3);
            pi_k    = zeros(n_landmarks+1, 1);
            
            % create temporary new landmark at observed position
            temp_mark = [stateMeanBar(1) + z(1,k) * cos(z(2,k) + stateMeanBar(3));
                         stateMeanBar(2) + z(1,k) * sin(z(2,k) + stateMeanBar(3));
                         0]; % z(3,k)]; % not sure whether to include landmark signature
            stateMeanTemp = [stateMeanBar;
                             temp_mark];
            stateCovTemp = [stateCovBar zeros(size(stateCovBar,2),1);
                            zeros(1,size(stateCovBar,2) + 1)];
            % initialize covariance for new landmark
            stateCovTemp(n_landmarks + 4, n_landmarks + 4) = 0.65;
            
            % loop over all landmarks  (including the temp landmark)and 
            % compute likelihood of correspondence with the new landmark
            % NOTE: could improve by caching predicted observations when
            %       more than 1 observation occurs at the same timestep
            for j = 1:n_landmarks+1

                delta = [mu_jx - stateMeanTemp(1);
                         mu_jy - stateMeanTemp(2)];
                       
                q = delta' * delta;
                r = sqrt(q);
                
                predZ(j) = [r;
                            conBear(atan2(delta(2), delta(1)) - stateMeanTemp(3));
                            0];
                F_xk = [eye(3)     zeros(3*j-3) zeros(3,3) zeros(3*n_landmarks+1 - 3*j);
                        zeros(3,3) zeros(3*j-3) eye(3)     zeros(3*n_landmarks+1 - 3*j)];
            end
        end
    end
end