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
        
% reshapes stateMean to 1 x 3*n_landmarks+3
F_x = [eye(3) zeros(3, 3 * n_landmarks)];

stateMean = F_x' * stateMean;

% initialize pose covariance with small, nonzero values  
stateCov = zeros(3 * n_landmarks + 3, 3 * n_landmarks + 3);
stateCov(1:3,1:3) = 0.001;

% initialize landmark covariance with known uncertainties
% on landmark positions
for i = 1:n_landmarks * 3 + 3
    %{
    stateCov(i * 3 + 1, i * 3 + 1) = Landmark_Groundtruth(i, 4);
    stateCov(i * 3 + 2, i * 3 + 2) = Landmark_Groundtruth(i, 5);
    %}
    stateCov(i,i) = 0.65;
end

    
measurementIndex = 1;

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
    stateMeanBar = stateMean + F_x' * poseUpdate;
    
    % calculate movement jacobian
    g_t = [0 0 trans * -sin(theta + halfRot);
                         0 0 trans * cos(theta + halfRot);
                         0 0 0];
    G_t = eye(3 * n_landmarks + 3, 3 * n_landmarks + 3) + F_x' * g_t * F_x;
     
    % calculate motion covariance in control space
    M_t = [(alphas(1) * abs(u_t(1)) + alphas(2) * abs(u_t(2)))^2 0;
           0 (alphas(3) * abs(u_t(1)) + alphas(4) * abs(u_t(2)))^2];
       
    % calculate Jacobian to transform motion covariance to state space
    V_t = [cos(theta + halfRot) -0.5 * sin(theta + halfRot);
           sin(theta + halfRot) 0.5 * cos(theta + halfRot);
           0 1];
    
    % update state covariance
    % not totally sure I'm calculating R_t correctly
    R_t = V_t * M_t * V_t';
    stateCovBar = (G_t * stateCov * G_t') + (F_x' * R_t * F_x);
    
    
    % build vector of features observed at current time
    z = zeros(3,1);
    while (Robots{robot_num}.M(measurementIndex, 1) - t < .005) && (measurementIndex < size(Robots{robot_num}.M,1))
        barcode = Robots{robot_num}.M(measurementIndex,2);
        landmarkID = 0;
        if (codeDict.isKey(barcode))
            landmarkID = codeDict(barcode);
        else
            disp('key not found');
        end
        if landmarkID > 5 && landmarkID < 21
            range = Robots{robot_num}.M(measurementIndex, 3);
            bearing = Robots{robot_num}.M(measurementIndex, 4);
            if uint8(z(3)) == 0
                z = [range;
                     bearing;
                     landmarkID - 5];
            else
                newZ = [range;
                        bearing;
                        landmarkID - 5];
                z = [z newZ];
            end
        end
        measurementIndex = measurementIndex + 1;
    end
    
    % if features are observed
    % loop over all features and compute Kalman gain
    if z(3,1) > 1
        
        S = zeros(size(z,2),3,3);
        zHat = zeros(3, size(z,2));
    
        for k = 1:size(z, 2) % loop over every observed landmark
            j = z(3,k);
            
            % if the landmark has never been seen before
            % add it to the state vector
            if stateMeanBar(3 + j) == 0
                landmark_pos = [z(1,k) * cos(z(2,k) + stateMeanBar(3));
                                z(1,k) * sin(z(2,k) + stateMeanBar(3));
                                0];
                stateMeanBar(3*j+1:3*j+3) = [stateMeanBar(1) + landmark_pos(1);
                                             stateMeanBar(2) + landmark_pos(2);
                                             0];
            end

            % compute predicted range and bearing
            delta = [stateMeanBar(3*j+1) - stateMeanBar(1);
                     stateMeanBar(3*j+2) - stateMeanBar(2)];
            q = delta' * delta;
            r = sqrt(q); % predicted range to landmark
            
            % predicted bearing to landmark
            pred_bear = conBear(atan2(delta(2), delta(1)) - stateMeanBar(3));
            
            zHat(:,k) = [r;
                         pred_bear;
                         j];
                     
            h_t = [-r*delta(1) -r*delta(2)  0   r*delta(1) r*delta(2) 0;
                   delta(2)    -delta(1)    -q  -delta(2)  delta(1)   0;
                   0           0            0   0          0          q];
                     
            F_1 = [eye(3,3); zeros(3,3)];
            F_2 = [zeros(3,3); eye(3,3)];
            F_xj = [F_1 zeros(6,3*j-3) F_2 zeros(6,3*n_landmarks - 3*j)];
            
            H_t = (1/q) * h_t * F_xj;
            
            % compute Kalman gain
            K = stateCovBar * H_t' * inv(H_t * stateCovBar * H_t' + Q_t);
            
            % incorporate new measurement into state mean and covariance
            stateMeanBar = stateMeanBar + K * (z(:,k) - zHat(:,k));
            stateCovBar = (eye(48) - (K * H_t)) * stateCovBar;
        end
    end
    
    % update state mean and covariance
    stateMean = stateMeanBar;
    stateCov = stateCovBar;
    
    % add new pose mean to estimated poses
    Robots{robot_num}.Est(i,:) = [t stateMean(1) stateMean(2) stateMean(3)];
end

% calculate land_loss: sum of squares of error in final landmark position
land_loss = 0;
for i = 1:n_landmarks
    x_diff = stateMean(3 * i + 1) - Landmark_Groundtruth(i, 2);
    y_diff = stateMean(3 * i + 2) - Landmark_Groundtruth(i, 3);
    sq_err = x_diff^2 + y_diff^2;
    land_loss = land_loss + sq_err;
end
%disp(land_loss);

p_loss = path_loss(Robots, robot_num, start);
disp(p_loss);
animateMRCLAMdataSet(Robots, Landmark_Groundtruth, timesteps, deltaT);