addpath ../common;

deltaT = .02;
alphas = [.1 .01 .18 .08 0 0]; % motion model noise parameters

% measurement model noise parameters
Q_t = [11.8 0;
       0    .18];
   
n_robots = 1;
robot_num = 1;
n_landmarks = 15;
[Barcodes, Landmark_Groundtruth, Robots] = loadMRCLAMdataSet(n_robots);
[Robots, timesteps] = sampleMRCLAMdataSet(Robots, deltaT);

% add pose estimate matrix to Robots
Robots{robot_num}.Est = zeros(size(Robots{robot_num}.G,1), 4);

start = 600;
t = Robots{robot_num}.G(start, 1);

% initialize state mean with the groundtruth
% pose at timestep zero
stateMean = [Robots{robot_num}.G(start,2);
            Robots{robot_num}.G(start,3);
            Robots{robot_num}.G(start,4)];
        
% reshapes stateMean to 1x(2*n_landmarks+3)
% and remaps motion updates to new state vector
F_x = [eye(3) zeros(3, 2 * n_landmarks)];

stateMean = F_x' * stateMean;

% initialize diagonal of pose covariance with small, nonzero values
% since we have a good estimate of the starting pose
stateCov = zeros(2 * n_landmarks + 3, 2 * n_landmarks + 3);
stateCov(1:3,1:3) = 0.001;

% initialize landmark variances to large values since we have
% no prior information on the locations of the landmarks
for i = 4:n_landmarks * 2 + 3
    stateCov(i,i) = 10;
end
   
measurementIndex = 1;

% set up map between barcodes and landmark IDs
% barcode values are NOT the same as landmark IDs
codeDict = containers.Map(Barcodes(:,2),Barcodes(:,1));

% loop through all odometry and measurement samples
% updating the robot's pose estimate with each step
% reference table 10.1 in Probabilistic Robotics
for i = start:size(Robots{robot_num}.G, 1)
    
    % update time
    t = Robots{robot_num}.G(i, 1);
    
    % get control inputs at current timestep
    u_t = Robots{robot_num}.O(i, 2:end);
    % update robot bearing to last timestep's estimate
    theta = stateMean(3, 1);
    % calculate rotation and translation over last timestep
    rot = deltaT * u_t(2);
    halfRot = rot / 2;
    trans = u_t(1) * deltaT;
    
    % calculate pose update in world frame
    poseUpdate = [trans * cos(theta + halfRot);
                  trans * sin(theta + halfRot);
                  rot];
              
    % updated state mean and constrain bearing
    % to +/- pi
    stateMeanBar = stateMean + F_x' * poseUpdate;
    stateMeanBar(3) = conBear(stateMeanBar(3));
    
    % calculate movement jacobian
    g_t = [0 0 trans * -sin(theta + halfRot);
           0 0 trans * cos(theta + halfRot);
           0 0 0];
    G_t = eye(2 * n_landmarks + 3) + F_x' * g_t * F_x;
     
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
    
    % get measurements for current timestep
    [z, measurementIndex] = getObservations(Robots, robot_num, t, measurementIndex, codeDict);

    % if measurements are available
    if z(3,1) > 1
        S = zeros(size(z,2),3,3);
        zHat = zeros(2, size(z,2));
        % loop over every measurement and use to correct state estimate
        for k = 1:size(z, 2) 
            j = z(3,k);
            
            % if the landmark has never been seen before
            % add it to the state vector
            if stateMeanBar(j*2 + 2) == 0
                landmark_pos = [z(1,k) * cos(z(2,k) + stateMeanBar(3));
                                z(1,k) * sin(z(2,k) + stateMeanBar(3))];
                stateMeanBar(2*j+2:2*j+3) = [stateMeanBar(1) + landmark_pos(1);
                                             stateMeanBar(2) + landmark_pos(2)];
                continue;
            end

            % compute predicted range and bearing to the landmark
            delta = [stateMeanBar(2*j+2) - stateMeanBar(1);
                     stateMeanBar(2*j+3) - stateMeanBar(2)];
            q = delta' * delta;
            r = sqrt(q); % predicted range to landmark
            
            % predicted bearing to landmark
            pred_bear = conBear(atan2(delta(2), delta(1)) - stateMeanBar(3));
            % create predicted measurement
            zHat(:,k) = [r;
                         pred_bear];
            % calculate measurement Jacobian    
            h_t = [-delta(1)/r -delta(2)/r  0   delta(1)/r delta(2)/r;
                   delta(2)/q    -delta(1)/q    -1  -delta(2)/q  delta(1)/q];
                     
            F_1 = [eye(3); zeros(2,3)];
            F_2 = [zeros(3,2); eye(2)];
            F_xj = [F_1 zeros(5,2*j-2) F_2 zeros(5,2*n_landmarks - 2*j)];
            
            H_t = h_t * F_xj;
            
            % compute Kalman gain
            K = stateCovBar * H_t' / (H_t * stateCovBar * H_t' + Q_t);
            
            % incorporate new measurement into state mean and covariance
            stateMeanBar = stateMeanBar + K * (z(1:2,k) - zHat(:,k));
            stateCovBar = (eye(2*n_landmarks+3) - (K * H_t)) * stateCovBar;
        end
    end
    
    % update state mean and covariance
    stateMean = stateMeanBar;
    stateMean(3) = conBear(stateMean(3));
    stateCov = stateCovBar;
    
    % add new pose mean to estimated poses
    Robots{robot_num}.Est(i,:) = [t stateMean(1) stateMean(2) stateMean(3)];
end

% calculate land_loss: sum of squares of error in final landmark position
land_loss = 0;
for i = 1:n_landmarks
    x_diff = stateMean(2 * i + 1) - Landmark_Groundtruth(i, 2);
    y_diff = stateMean(2 * i + 2) - Landmark_Groundtruth(i, 3);
    sq_err = x_diff^2 + y_diff^2;
    land_loss = land_loss + sq_err;
end
disp(land_loss);

p_loss = path_loss(Robots, robot_num, start);
disp(p_loss);
animateMRCLAMdataSet(Robots, Landmark_Groundtruth, timesteps, deltaT);