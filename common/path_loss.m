% computes euclidean loss between robot's estimated path and ground truth
% ignores bearing error
function loss = path_loss(Robots, robot_num, start)
    loss = 0;
    for i = start:size(Robots{robot_num}.G,1)
        x_diff = Robots{robot_num}.G(i,2) - Robots{robot_num}.Est(i,2);
        y_diff = Robots{robot_num}.G(i,3) - Robots{robot_num}.Est(i,3);
        err = x_diff^2 + y_diff^2;
        loss = loss + err;
    end
end