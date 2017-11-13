addpath ../common;

deltaT = .02;
alphas = [.5 .5 1 1 .5 .5]; % need to figure out how to set these

% also don't know how to calculate the measurement noise std_dev
sigma_range = 300;
sigma_bearing = 300;
sigma_id = 150;

Q_t = [sigma_range^2 0 0;
       0 sigma_bearing^2 0;
       0 0 sigma_id^2];

measurement_prob = 0;
n_robots = 1;
robot_num = 1;
n_landmarks = 15;
[Barcodes, Landmark_Groundtruth, Robots] = loadMRCLAMdataSet(n_robots);
[Robots, timesteps] = sampleMRCLAMdataSet(Robots, deltaT);

% add pose estimate matrix to Robots
Robots{robot_num}.Est = zeros(size(Robots{robot_num}.G,1), 4);

% initialize time, and pose estimate
start = 600; %15000;
t = Robots{robot_num}.G(start, 1);
% need mean and covariance for the initial pose estimate
poseMean = [Robots{robot_num}.G(start,2);
            Robots{robot_num}.G(start,3);
            Robots{robot_num}.G(start,4)];
poseCov = [0.01 0.01 0.01;
           0.01 0.01 0.01;
           0.01 0.01 0.01];
       
measurementIndex = 2;
angularCorrect = .01;

% set up map between barcodes and landmark IDs
codeDict = containers.Map(Barcodes(:,2),Barcodes(:,1));
%{
for i = 1:size(Barcodes,1)
    codeDict(Barcodes(i,2)) = Barcodes(i,1);
end
%}
while (Robots{robot_num}.M(measurementIndex, 1) < t - .05)
        measurementIndex = measurementIndex + 1;
end

% loop through all odometry and measurement samples
% updating the robot's pose estimate with each step
% reference table 7.2 in Probabilistic Robotics
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
    
    if (mod(uint32(t*100), 100) == 0)
        1;
    end
    
    % build vector of features observed at current time
    z = zeros(3,1);
    %
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
                     0];
            else
                newZ = [range;
                        bearing;
                        0];
                z = [z newZ];
            end
        end
        measurementIndex = measurementIndex + 1;
    end
    
    % if features are observed
    % loop over all features and compute Kalman gain
    if z(1,1) > 0.1
        for k = 1:size(z, 2) % loop over every observed landmark
            
            % loop over all landmarks and compute MLE correspondence
            predZ = zeros(n_landmarks, 1, 3);
            predS = zeros(n_landmarks, 3, 3);
            predH = zeros(n_landmarks, 3, 3);
            maxJ = 0;
            landmarkIndex = 0;
            for j = 1:n_landmarks
                xDist = Landmark_Groundtruth(j, 2) - poseMeanBar(1);
                yDist = Landmark_Groundtruth(j, 3) - poseMeanBar(2);
                q = xDist^2 + yDist^2;
                predZ(j,:,:) = [sqrt(q);
                            conBear(atan2(yDist, xDist) - poseMeanBar(3));
                            0];
                predH(j,:,:) = [-xDist/sqrt(q) -yDist/sqrt(q) 0;
                                yDist/q        -xDist/q      -1;
                                0              0              0];
                predS(j,:,:) = predH(j) * poseCovBar * predH(j)' + Q_t;
                thisJ = det(2 * pi * squeeze(predS(j,:,:)))^(-0.5) * ...
                        exp(-0.5 * (z(:,k) - squeeze(predZ(j,:,:)))' ...
                        * inv(squeeze(predS(j,:,:))) ...
                        * (z(:,k) - squeeze(predZ(j,:,:))));
                if imag(thisJ) ~= 0
                    thisJ = 0;
                end
                if thisJ > maxJ
                    maxJ = thisJ;
                    landmarkIndex = j;
                elseif thisJ < 0
                    disp("j less than 0");
                    disp(thisJ);
                elseif thisJ == 0
                    disp("j equals 0");
                end
            end
            
            % compute Kalman gain
            K = poseCovBar * squeeze(predH(landmarkIndex,:,:))' ...
                * inv(squeeze(predS(landmarkIndex,:,:)));
            %{
            xAct = m(1)-Robots{robot_num}.G(i,2);
            yAct = m(2)-Robots{robot_num}.G(i,3);
            zAct = [sqrt(xAct^2 + yAct^2);
                    atan2(yAct, xAct) - Robots{robot_num}.G(i,4);
                    j];
            %}
            % update pose mean and covariance estimates
            poseMeanBar = poseMeanBar + K * (z(:,k) - squeeze(predZ(landmarkIndex,:,:)));
            poseCovBar = (eye(3) - (K * squeeze(predH(landmarkIndex,:,:)))) * poseCovBar;
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
    Robots{robot_num}.Est(i,:) = [t poseMean(1) poseMean(2) poseMean(3)];
    % calculate error between mean pose and groundtruth
    %{
    groundtruth = [Robots{robot_num}.G(i, 2); 
                   Robots{robot_num}.G(i, 3); 
                   Robots{robot_num}.G(i, 4)];
    error = groundtruth - poseMean;
    %}
   
end

animateMRCLAMdataSet(Robots, Landmark_Groundtruth, timesteps, deltaT);