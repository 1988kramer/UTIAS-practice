function landmarks_in_FOV = in_FOV(predZ, maxBear, minRange, maxRange)
    n_landmarks = size(predZ, 2);
    landmarks_in_FOV = zeros(n_landmarks, 1);
    for i = 1:n_landmarks
        if predZ(1,i) < maxRange ...
                && predZ(1,i) > minRange ...
                && (abs(predZ(2,i)) < maxBear)
            landmarks_in_FOV(i) = 1;
        end
    end
end