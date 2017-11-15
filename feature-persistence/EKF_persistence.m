addpath ../common;

deltaT = .02;
alphas = [.2 .03 .09 .08 0 0]; % motion model noise parameters

% measurement model noise parameters
sigma_range = .43;
sigma_bearing = .6;
sigma_id = 1;

Q_t = [sigma_range^2 0 0;
       0 sigma_bearing^2 0;
       0 0 sigma_id^2];
   
n_robots = 1;
robot_num = 1;
n_landmarks = 15;
   
pM = 0.6;      % probability of missed measurement
pF = 0.3;      % probability of false positive measurement
lambda = 1e-2; % prior survival function = exp(-lambda * t)
pV = 1e-10;    % threshold probability for landmark persistence

n_incorrect = 0;       % number of incorrect landmark associations
misses = zeros(n_landmarks, 1);   % predicted measurements for which a landmark was 
                       % not detected when it is in the sensor's FOV
expected = zeros(3,1); % predicted measurements which were in 
                       % the sensor's FOV
[Barcodes, Landmark_Groundtruth, Robots] = loadMRCLAMdataSet(n_robots);
[Robots, timesteps] = sampleMRCLAMdataSet(Robots, deltaT);

% persist(:,1) = posterior persistence probability
% persist(:,2) = partial evidence
% persist(:,3) = likelihood
% persist(:,4) = evidence
% persist(:,5) = persistence prior
% persist(:,6) = time landmark was first observed
persist = zeros(n_landmarks,6); % may need to change how this is initialized
persist(:,3) = 1;
persist(:,4) = 1;

% add pose estimate matrix to Robots
Robots{robot_num}.Est = zeros(size(Robots{robot_num}.G,1), 4);

start = 1;
t = Robots{robot_num}.G(start, 1);
tmax = max(Robots{robot_num}.G(:,1));

% need mean and covariance for the initial pose estimate
poseMean = [Robots{robot_num}.G(start,2);
            Robots{robot_num}.G(start,3);
            Robots{robot_num}.G(start,4)];
poseCov = [0.01 0.01 0.01;
           0.01 0.01 0.01;
           0.01 0.01 0.01];
       
measurementIndex = 1;

% set up map between barcodes and landmark IDs
codeDict = containers.Map(Barcodes(:,2),Barcodes(:,1));

% randomly decide times when several landmarks will disappear
times_of_death = kill_landmarks(n_landmarks, 2, t, 30, deltaT);

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
    
    % stopping point for debugging
    if (mod(uint32(t*100), 100) == 0)
        1;
    end
    
    % build vector of features observed at current time
    z = zeros(3,1);
    landmarkIDs = 0;
    while (Robots{robot_num}.M(measurementIndex, 1) - t < .005) && (measurementIndex < size(Robots{robot_num}.M,1))
        barcode = Robots{robot_num}.M(measurementIndex,2);
        landmarkID = 0;
        if (codeDict.isKey(barcode))
            landmarkID = codeDict(barcode);
        else
            disp('key not found');
        end
        % check if landmark is valid
        if landmarkID > 5 && landmarkID < 21
            TOD = times_of_death(landmarkID - 5);
            % check if landmark has died
            if (TOD == 0 || TOD > t)
                range = Robots{robot_num}.M(measurementIndex, 3);
                bearing = Robots{robot_num}.M(measurementIndex, 4);
                if z(3,1) == 0
                    z = [range;
                         bearing;
                         0];
                    landmarkIDs = landmarkID;
                else
                    newZ = [range;
                            bearing;
                            0];
                    z = [z newZ];
                    newID = landmarkID;
                    landmarkIDs = [landmarkIDs newID];
                end
            %{
            else
                misses(landmarkID - 5) = misses(landmarkID - 5) + 1;
            %}
            end
        end
        measurementIndex = measurementIndex + 1;
    end
    
    % if features are observed
    % loop over all features and compute Kalman gain
    if z(1,1) > 0.1
        predZ = zeros(3, n_landmarks);
        for k = 1:size(z, 2) % loop over every observed landmark
            
            % loop over all landmarks and compute MLE correspondence
            predS = zeros(n_landmarks, 3, 3);
            predH = zeros(n_landmarks, 3, 3);
            maxJ = 0;
            landmarkIndex = 0;
            for j = 1:n_landmarks
                xDist = Landmark_Groundtruth(j, 2) - poseMeanBar(1);
                yDist = Landmark_Groundtruth(j, 3) - poseMeanBar(2);
                q = xDist^2 + yDist^2;
                predZ(:,j) = [sqrt(q);
                            conBear(atan2(yDist, xDist) - poseMeanBar(3));
                            0];
                predH(j,:,:) = [-xDist/sqrt(q) -yDist/sqrt(q) 0;
                                yDist/q        -xDist/q      -1;
                                0              0              0];
                predS(j,:,:) = squeeze(predH(j,:,:)) * poseCovBar * ...
                               squeeze(predH(j,:,:))' + Q_t;
                thisJ = det(2 * pi * squeeze(predS(j,:,:)))^(-0.5) * ...
                        exp(-0.5 * (z(:,k) - predZ(:,j))' ...
                        * inv(squeeze(predS(j,:,:))) ...
                        * (z(:,k) - predZ(:,j)));

                if imag(thisJ) ~= 0
                    thisJ = 0;
                end
                if thisJ > maxJ && persist(j,6) >= 0
                    maxJ = thisJ;
                    landmarkIndex = j;
                elseif thisJ < 0
                    disp("j less than 0");
                    disp(thisJ);
                elseif thisJ == 0
                    disp("j equals 0");
                end
            end
            z(3,k) = landmarkIndex;
            
            % increment missed count if measurement is associated
            % to the wrong landmark
            if landmarkIDs(k) ~= landmarkIndex + 5
                n_incorrect = n_incorrect + 1;
            end
            
            % if landmark persistence probability is above threshold
            % update pose mean and covariance with measurement
            % else ignore the measurement
            % NOTE persistence probabilities should be updated before this
            % compute Kalman gain
            K = poseCovBar * squeeze(predH(landmarkIndex,:,:))' ...
                * inv(squeeze(predS(landmarkIndex,:,:)));

            % update pose mean and covariance estimates
            poseMeanBar = poseMeanBar + K * (z(:,k) - predZ(:,landmarkIndex));
            poseCovBar = (eye(3) - (K * squeeze(predH(landmarkIndex,:,:)))) * poseCovBar;
        end
        % update landmark persistence probabilities
        marks_in_FOV = in_FOV(predZ, 0.3, 1, 4);
        [persist, misses] = persistence_filter(pM, pF, lambda, marks_in_FOV, ...
                                         persist, z, t, pV, misses);
    end
    
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

%animateMRCLAMdataSet(Robots, Landmark_Groundtruth, timesteps, deltaT);
