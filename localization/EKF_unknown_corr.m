addpath ../common;

deltaT = .02;
alphas = [.2 .03 .09 .08 0 0]; % motion model noise parameters

% measurement model noise parameters
sigma_range = .43;
sigma_bearing = .6;

Q_t = [sigma_range^2 0;
       0 sigma_bearing^2];

n_robots = 1;
robot_num = 1;
n_landmarks = 15;
[Barcodes, Landmark_Groundtruth, Robots] = loadMRCLAMdataSet(n_robots);
[Robots, timesteps] = sampleMRCLAMdataSet(Robots, deltaT);

% add pose estimate matrix to Robots
Robots{robot_num}.Est = zeros(size(Robots{robot_num}.G,1), 4);

start = 1; 
t = Robots{robot_num}.G(start, 1); 
% need mean and covariance for the initial pose estimate
poseMean = [Robots{robot_num}.G(start,2);
            Robots{robot_num}.G(start,3);
            Robots{robot_num}.G(start,4)];
poseCov = 0.01*eye(3);
       
measurementIndex = 1;

% set up map between barcodes and landmark IDs
codeDict = containers.Map(Barcodes(:,2),Barcodes(:,1));

while (Robots{robot_num}.M(measurementIndex, 1) < t - .05)
        measurementIndex = measurementIndex + 1;
end

% loop through all odometry and measurement samples
% updating the robot's pose estimate with each step
% reference table 7.3 in Probabilistic Robotics
for i = start:size(Robots{robot_num}.G, 1)
    theta = poseMean(3, 1);
    % update time
    t = Robots{robot_num}.G(i, 1);
    % update movement vector
    u_t = [Robots{robot_num}.O(i, 2); Robots{robot_num}.O(i, 3)];

    rot = deltaT * u_t(2);
    halfRot = rot / 2;
    trans = u_t(1) * deltaT;
    
    % calculate the movement Jacobian
    G_t = [1 0 trans * -sin(theta + halfRot);
           0 1 trans * cos(theta + halfRot);
           0 0 1];
    % calculate motion covariance in control space
    M_t = [(alphas(1) * abs(u_t(1)) + alphas(2) * abs(u_t(2)))^2 0;
           0 (alphas(3) * abs(u_t(1)) + alphas(4) * abs(u_t(2)))^2];
    % calculate Jacobian to transform motion covariance to state space
    V_t = [cos(theta + halfRot) -0.5 * sin(theta + halfRot);
           sin(theta + halfRot) 0.5 * cos(theta + halfRot);
           0 1];
    
    % calculate pose update from odometry
    poseUpdate = [trans * cos(theta + halfRot);
                  trans * sin(theta + halfRot);
                  rot];
    
    poseMeanBar = poseMean + poseUpdate;
    poseCovBar = G_t * poseCov * G_t' + V_t * M_t * V_t';
    
    % get measurements
    [z, measurementIndex] = getObservations(Robots, robot_num, t, measurementIndex, codeDict);
    
    % remove landmark ID from measurement because
    % we're trying to live without that
    z = z(1:2,:);
    
    % if features are observed
    % loop over all features and compute Kalman gain
    if z(1,1) > 0.1
        for k = 1:size(z, 2) % loop over every observed landmark
            
            % loop over all landmarks and compute MLE correspondence
            predZ = zeros(n_landmarks, 1, 2); % predicted measurements
            predS = zeros(n_landmarks, 2, 2); % predicted measurement covariances
            predH = zeros(n_landmarks, 2, 3); % predicted measurement Jacobians
            maxJ = 0; 
            landmarkIndex = 0;
            for j = 1:n_landmarks
                xDist = Landmark_Groundtruth(j, 2) - poseMeanBar(1);
                yDist = Landmark_Groundtruth(j, 3) - poseMeanBar(2);
                q = xDist^2 + yDist^2;
                
                % calculate predicted measurement
                predZ(j,:,:) = [sqrt(q);
                            conBear(atan2(yDist, xDist) - poseMeanBar(3))];
                        
                % calculate predicted measurement Jacobian
                predH(j,:,:) = [-xDist/sqrt(q) -yDist/sqrt(q) 0;
                                yDist/q        -xDist/q      -1];
                            
                % calculate predicted measurement covariance
                predS(j,:,:) = squeeze(predH(j,:,:)) * poseCovBar ...
                               * squeeze(predH(j,:,:))' + Q_t;
                
                % calculate probability of measurement correspondence          
                thisJ = det(2 * pi * squeeze(predS(j,:,:)))^(-0.5) * ...
                        exp(-0.5 * (z(:,k) - squeeze(predZ(j,:,:)))' ...
                        * inv(squeeze(predS(j,:,:))) ...
                        * (z(:,k) - squeeze(predZ(j,:,:))));
                
                % update correspondence if the probability is
                % higher than the previous maximum
                if thisJ > maxJ
                    maxJ = thisJ;
                    landmarkIndex = j;
                end
            end

            % compute Kalman gain
            K = poseCovBar * squeeze(predH(landmarkIndex,:,:))' ...
                * inv(squeeze(predS(landmarkIndex,:,:)));

            % update pose mean and covariance estimates
            poseMeanBar = poseMeanBar + K * (z(:,k) - squeeze(predZ(landmarkIndex,:,:)));
            poseCovBar = (eye(3) - (K * squeeze(predH(landmarkIndex,:,:)))) * poseCovBar;
        end
    end
    
    % update pose mean and covariance
    poseMean = poseMeanBar;
    poseMean(3) = conBear(poseMean(3));
    poseCov = poseCovBar;
    
    % add pose mean to estimated position vector
    Robots{robot_num}.Est(i,:) = [t poseMean(1) poseMean(2) poseMean(3)];
    
    % calculate error between mean pose and groundtruth 
    % for testing only
    %{
    groundtruth = [Robots{robot_num}.G(i, 2); 
                   Robots{robot_num}.G(i, 3); 
                   Robots{robot_num}.G(i, 4)];
    error = groundtruth - poseMean;
    %}
   
end
animateMRCLAMdataSet(Robots, Landmark_Groundtruth, timesteps, deltaT);