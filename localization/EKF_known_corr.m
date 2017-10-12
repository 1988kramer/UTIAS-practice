addpath ../common;

deltaT = .02;
alphas = [1 1 10 10 1 1]; % need to figure out how to set these

% also don't know how to calculate the measurement noise std_dev
sigma_range = 75;
sigma_bearing = 125;
sigma_id = 50;

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
start = 15000;
t = Robots{1}.G(start, 1);
% need mean and covariance for the initial pose estimate
poseMean = [Robots{1}.G(start,2);
            Robots{1}.G(start,3);
            Robots{1}.G(start,4)];
poseCov = [0.01 0.01 0.01;
           0.01 0.01 0.01;
           0.01 0.01 0.01];
       
measurementIndex = 2;
angularCorrect = .01;

while (Robots{1}.M(measurementIndex, 1) < t - .05)
        measurementIndex = measurementIndex + 1;
end

% loop through all odometry and measurement samples
% updating the robot's pose estimate with each step
% reference table 7.2 in Probabilistic Robotics
for i = start:size(Robots{1}.G, 1)
    theta = poseMean(3, 1);
    % update time
    t = Robots{1}.G(i, 1);
    % update movement vector
    u_t = [Robots{1}.O(i, 2); Robots{1}.O(i, 3)];
    
    if u_t(2) == 0
        u_t(2) = angularCorrect;
        angularCorrect = angularCorrect * -1;
    end
    
    % calculate a bunch of commonly used terms
    sin1 = sin(theta + (u_t(2) * deltaT));
    cos1 = cos(theta + (u_t(2) * deltaT));
    plusSin = (sin(theta) + sin1);
    plusCos = (cos(theta) + cos1);
    minusSin = (sin(theta) - sin1);
    minusCos = (cos(theta) - cos1);
    
    % calculate the movement Jacobian
    G_t = [1 0 (u_t(1) / u_t(2)) * (-1 * minusCos);
           0 1 (u_t(1) / u_t(2)) * (-1 * minusSin);
           0 0 1];
    % calculate motion covariance in control space
    M_t = [(alphas(1) * abs(u_t(1)) + alphas(2) * abs(u_t(2)))^2 0;
           0 (alphas(3) * abs(u_t(1)) + alphas(4) * abs(u_t(2)))^2];
    % calculate Jacobian to transform motion covariance to state space
    V_t = [(-1*plusSin/u_t(2)) (((u_t(1)*minusSin)/(u_t(2)^2))+(u_t(1)*cos1*deltaT)/u_t(2));
           (minusCos/u_t(2)) (((-1*u_t(1)*minusCos)/(u_t(2)^2))+(u_t(1)*sin1*deltaT)/u_t(2));
           0 deltaT];
    %{
    % calculate pose update from odometry
    poseUpdate = [u_t(1) * cos(theta) * deltaT;
                  u_t(1) * sin(theta) * deltaT;
                  u_t(2) * deltaT];
    %}
    poseUpdate = [-1 * (u_t(1) / u_t(2)) * minusSin;
                  (u_t(1) / u_t(2)) * minusCos;
                  u_t(2) * deltaT];
    poseMeanBar = poseMean + poseUpdate;
    poseCovBar = G_t * poseCov * G_t' + V_t * M_t * V_t';
    
    if (mod(uint32(t*100), 100) == 0)
        1;
    end
    
    % build vector of features observed at current time
    z = zeros(3,1);
    %
    while (Robots{1}.M(measurementIndex, 1) - t < .005) && (measurementIndex < size(Robots{1}.M,1))
        landmarkID = Robots{1}.M(measurementIndex,2);
        if landmarkID > 5 && landmarkID < 21
            range = Robots{1}.M(measurementIndex, 3);
            bearing = Robots{1}.M(measurementIndex, 4);
            id = landmarkID;
            if uint8(z(3)) == 0
                z = [range;
                     bearing;
                     id];
            else
                newZ = [range;
                        bearing;
                        id];
                z = [z newZ];
            end
        end
        measurementIndex = measurementIndex + 1;
    end
    S = zeros(size(z,2),3,3);
    zHat = zeros(3, size(z,2));
    
    % if features are observed
    % loop over all features and compute Kalman gain
    if z(3,1) > 1
        for k = 1:size(z, 2) % loop over every observed landmark
            j = z(3,k);

            % get coordinates of observed landmark
            m = Landmark_Groundtruth(j - 5, 2:3);

            % compute predicted range and bearing
            xDist = m(1) - poseMeanBar(1);
            yDist = m(2) - poseMeanBar(2);
            q = xDist^2 + yDist^2;
            
            pred_bear = conBear(atan2(yDist, xDist) - poseMeanBar(3));
            
            zHat(:,k) = [sqrt(q);
                         pred_bear;
                         j];

            % calculate Jacobian of h (line 13)
            H = [(-1 * (xDist / sqrt(q))) (-1 * (yDist / sqrt(q))) 0;
                 (yDist / q) (-1 * (xDist / q)) -1
                 0 0 0];
            S(k,:,:) = H * poseCovBar * H' + Q_t;

            % compute Kalman gain
            K = poseCov * H' * inv(squeeze(S(k,:,:)));
            %K = squeeze(S(k,:,:)) \ (poseCovBar * H'); % may be equivalent to above

            % update pose mean and covariance estimates
            poseMeanBar = poseMeanBar + K * (z(:,k) - zHat(:,k));
            poseCovBar = (eye(3) - (K * H)) * poseCovBar;
        end
    end
    %
    % update pose mean and covariance
    poseMean = poseMeanBar;
    poseMean(3) = conBear(poseMean(3));
    poseCov = poseCovBar;
    %{
    % calculate measurement probability
    measurementProb = 1;
    if z(3,1) > 1
        for k = 1:size(z,2)
            detS = det(2 * pi * squeeze(S(k,:,:)))^(-0.5);
            expErr = exp(-0.5 * (z(:,k) - zHat(:,k))' * inv(squeeze(S(k,:,:))) * (z(:,k) - zHat(:,k)));
            %expErr = exp(-0.5 * (squeeze(S(k,:,:)) \ (z(:,k) - zHat(:,k))') * (z(:,k) - zHat(:,k)));
            measurementProb = measurementProb * detS * expErr;
        end
    end
    %}
    
    % add pose mean to estimated position vector
    Robots{1}.Est(i,:) = [t poseMean(1) poseMean(2) poseMean(3)];
    % calculate error between mean pose and groundtruth
    %{
    groundtruth = [Robots{1}.G(i, 2); 
                   Robots{1}.G(i, 3); 
                   Robots{1}.G(i, 4)];
    error = groundtruth - poseMean;
    %}
   
end

animateMRCLAMdataSet(Robots, Landmark_Groundtruth, timesteps, deltaT);