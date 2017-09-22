addpath ../common;

deltaT = .02;
alphas = [.01 .01 .01 .01 .01 .01]; % need to figure out how to set these

% also don't know how to calculate the measurement noise std_dev
sigma_range = .01;
sigma_bearing = .01;
sigma_id = 0;

Q_t = [sigma_range^2 0 0;
       0 sigma_bearing^2 0;
       0 0 sigma_id^2];

measurement_prob = 0;
n_robots = 1;
[Barcodes, Landmark_Groundtruth, Robots] = loadMRCLAMdataSet(n_robots);
[Robots, timesteps] = sampleMRCLAMdataSet(Robots, deltaT);

% add pose estimate matrix to Robots
Robots{1}.Est = zeros(size(Robots{1}.G,1), 4);

% initialize time, and pose estimate
t = 0;

% need mean and covariance for the initial pose estimate
poseMean = [Robots{1}.G(1,2);
            Robots{1}.G(1,3);
            Robots{1}.G(1,4)];
poseCov = [.001 .001 .001;
           .001 .001 .001;
           .001 .001 .001];
       
measurementIndex = 1;

% loop through all odometry and measurement samples
% updating the robot's pose estimate with each step
% reference table 7.2 in Probabilistic Robotics
for i = 1:size(Robots{1}.G, 1)
    theta = poseMean(3, 1);
    % update time
    t = Robots{1}.O(i, 1);
    % update movement vector
    u_t = [Robots{1}.O(i, 2); Robots{1}.O(i, 3)];
    
    % calculate a bunch of commonly used terms
    sin1 = sin(theta + u_t(2) * deltaT);
    cos1 = cos(theta + u_t(2) * deltaT);
    plusSin = (sin(theta) + sin1);
    plusCos = (cos(theta) + cos1);
    minusSin = (sin(theta) - sin1);
    minusCos = (cos(theta) - cos1);
    
    % calculate the movement Jacobian
    G_t = [1 0 (u_t(1) / u_t(2)) * ((-1 * cos(theta)) + plusCos);
           0 1 (u_t(1) / u_t(2)) * ((-1 * sin(theta)) + plusSin);
           0 0 1];
    % calculate motion covariance in control space
    M_t = [(alphas(1) * abs(u_t(1)) + alphas(2) * abs(u_t(2)))^2 0;
           0 (alphas(3) * abs(u_t(1)) + alphas(4) * abs(u_t(2)))^2];
    % calculate Jacobian to transform motion covariance to state space
    V_t = [(-1*plusSin/u_t(2)) (((u_t(1)*minusSin)/(u_t(2)^2))+(u_t(1)*cos1*deltaT)/u_t(2));
           (minusCos/u_t(2)) (((-1*minusCos)/(u_t(2)^2))+(u_t(1)*sin1*deltaT)/u_t(2));
           0 deltaT];
    
    % calculate pose update from odometry
    poseUpdate = [(u_t(1) / u_t(2)) * ((-1 * sin(theta)) + sin1);
                  (u_t(1) / u_t(2)) * (cos(theta) - cos1)
                  u_t(2) * deltaT];
    poseMeanBar = poseMean + poseUpdate;
    poseCovBar = G_t * poseCov * G_t' + V_t * M_t * V_t';
    
    
    
    % build vector of features observed at current time
    if Robots{1}.M(measurementIndex, 1) == t && Robots{i}.M(measurementIndex, 2) > 5
        range = Robots{1}.M(measurementIndex, 3);
        bearing = Robots{1}.M(measurementIndex, 4);
        id = Robots{1}.M(measurementIndex, 2);
        z = [range;
             bearing;
             id];
        measurementIndex = measurementIndex + 1;
    
        % record additional measurements if more than 
        % one happens at time t
        while Robots{1}.M(measurementIndex, 1) == t && Robots{i}.M(measurementIndex, 2) > 5
            range = Robots{1}.M(measurementIndex, 3);
            bearing = Robots{1}.M(measurementIndex, 4);
            id = Robots{1}.M(measurementIndex, 2);
            nextCol = [range;
                       bearing;
                       id];
            z = [z nextCol];
            measurementIndex = measurementIndex + 1;
        end
        S = zeros(size(z,2),3,3);
        zHat = zeros(3, size(z,2));

        % loop over all features and compute Kalman gain
        for k = 1:size(z, 2) % loop over every observed landmark
            j = z(3,k);

            % get coordinates of observed landmark
            m = Landmark_Groundtruth(j - 5, 2:3);

            % compute predicted range and bearing
            xDist = m(1) - poseMeanBar(1);
            yDist = m(2) - poseMeanBar(2);
            q = xDist^2 + yDist^2;
            zHat(:,k) = [sqrt(q);
                         atan2(yDist, xDist) - poseMean(3);
                         j];

            % calculate Jacobian of h (line 13)
            H = [(-1 * (xDist / sqrt(q))) (-1 * (yDist / sqrt(q))) 0;
                 (yDist / q) (-1 * (xDist / q)) -1
                 0 0 0];
            S(k,:,:) = H * poseCov * H' + Q_t;

            % compute Kalman gain
            % K = poseCov * H' * inv(S);
            K = S(k,:,:) \ (poseCov * H'); % may be equivalent to above

            % update pose mean and covariance estimates
            poseMeanBar = poseMeanBar + K * (z(:,k) - zHat(:,k));
            poseCovBar = (eye(3) - (K * H)) * poseCovBar;
        end
        
            % calculate measurement probability
        measurementProb = 1;
        for k = 1:size(z,2)
            detS = det(2 * pi * S(k,:,:))^(-0.5);
            % expErr = exp(-0.5 * (z(i,:,:) - zHat(i,:,:))' * inv(S(i,:,:)) * (z(i,:,:) - zHat(i,:,:)));
            expErr = exp(-0.5 * (S(k,:,:) \ (z(:,k) - zHat(:,k))') * (z(:,k) - zHat(:,k)));
            measurementProb = measurementProb * detS * expErr;
        end
    end
    
    % update pose mean and covariance
    poseMean = poseMeanBar;
    poseCov = poseCovBar;
    

    % add pose mean to estimated position vector
    Robots{1}.Est(i,:) = [t poseMean(1) poseMean(2) poseMean(3)];
    % calculate error between mean pose and groundtruth
    groundtruth = [Robots{1}.G(i, 2); 
                   Robots{1}.G(i, 3); 
                   Robots{1}.G(i, 4)];
    error = groundtruth - poseMean;
end

animateMRCLAMdataSet(Robots, Landmark_Groundtruth, timesteps, deltaT);