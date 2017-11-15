% returns an array of randomly selected times to kill the specified 
% number of randomly chosen landmarks
% NOTE: death of certain, less used landmarks can go unnoticed
%       should think of better way to select landmarks to kill
function times_of_death = kill_landmarks(n_landmarks, n_killed, t0, tmax, deltaT)
    times_of_death = zeros(n_landmarks, 1); % times each landmark is killed
    to_kill = randperm(n_landmarks, n_killed); % indices of landmarks to kill
    kill_times = randi([int64(t0 / deltaT), int64(tmax / deltaT)], 1, n_killed);
    kill_times = double(kill_times) .* deltaT;
    for i = 1:n_killed
        times_of_death(to_kill(i)) = kill_times(i);
    end
end