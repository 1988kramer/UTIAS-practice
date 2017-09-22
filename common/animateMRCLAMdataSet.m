% UTIAS Multi-Robot Cooperative Localization and Mapping Dataset
% produced by Keith Leung (keith.leung@robotics.utias.utoronto.ca) 2009
% Matlab script animateMRCLAMdataSet.m
% Description: This scripts creates an animation using ground truth data.
% Run this script after loadMRCLAMdataSet.m and sampleMRCLAMdataSet.m

function animateMRCLAMdataSet(Robots, Landmark_Groundtruth, timesteps, sample_time)

    % Options %
    start_timestep = 1;
    end_timestep = timesteps; 
    timesteps_per_frame = 50;
    pause_time_between_frames=0.01; %[s]
    draw_measurements = 0;
    % Options END %

    n_robots = size(Robots, 1);
    n_landmarks = length(Landmark_Groundtruth(:,1));

    % Plots and Figure Setup
    % Robot Colors
    colour(1,:) = [1 0 0];
    colour(2,:) = [0 0.75 0];
    colour(3,:) = [0 0 1];
    colour(4,:) = [1 0.50 0.25];
    colour(5,:) = [1 0.5 1];
    % Estimated Robot Colors
    colour(6,:) = [.67 0 0];
    colour(7,:) = [0 0.5 0];
    colour(8,:) = [0 0 .67];
    colour(9,:) = [.67 0.33 0.17];
    colour(10,:) = [.67 0.33 .67];
    for i=11:11+n_landmarks
        colour(i,:) = [0.3 0.3 0.3]; 
    end

    figHandle = figure('Name','Dataset Groundtruth','Renderer','OpenGL');
    set(gcf,'Position',[1300 1 630 950])
    plotHandles_robot_gt = zeros(n_robots,1);
    plotHandles_robot_est = zeros(n_robots,1);
    plotHandles_landmark_gt = zeros(n_landmarks,1);

    r_robot = 0.165;
    d_robot = 2*r_robot;
    r_landmark = 0.055;
    d_landmark = 2*r_landmark;
    
    % get initial positions for robot groundtruth
    for i = 1:n_robots  
        x=Robots{i}.G(1,2);
        y=Robots{i}.G(1,3);
        z=Robots{i}.G(1,4);
        x1 = d_robot*cos(z) + x;
        y1 = d_robot*sin(z) + y;
        p1 = x - r_robot;
        p2 = y - r_robot;
        plotHandles_robot_gt(i) = rectangle('Position',[p1,p2,d_robot,d_robot],'Curvature',[1,1],...
                  'FaceColor',colour(i,:),'LineWidth',1);
        line([x x1],[y y1],'Color','k');    
        n_measurements(i) = length(Robots{i}.M(:,1));
    end
    
    % get initial positions for robot pose estimates
    for i = 1:n_robots
        x=Robots{i}.Est(1,2);
        y=Robots{i}.Est(1,3);
        z=Robots{i}.Est(1,4);
        x1 = d_robot*cos(z) + x;
        y1 = d_robot*sin(z) + y;
        p1 = x - r_robot;
        p2 = y - r_robot;
        plotHandles_robot_est(i) = rectangle('Position',[p1,p2,d_robot,d_robot],'Curvature',[1,1],...
                  'FaceColor',colour(i + 5,:),'LineWidth',1);
        line([x x1],[y y1],'Color','k');    
    end
    
    % get positions for landmarks
    for i = 1:n_landmarks
        x=Landmark_Groundtruth(i, 2);
        y=Landmark_Groundtruth(i, 3);
        p1 = x - r_landmark;
        p2 = y - r_landmark;
        plotHandles_landmark_gt(i) = rectangle('Position',[p1,p2,d_landmark,d_landmark],'Curvature',[1,1],...
                  'FaceColor',colour(i+10,:),'LineWidth',1);
    end

    axis square;
    axis equal;
    axis([-2 6 -6 7]);
    set(gca,'XTick',(-10:2:10)');

    % Going throuhg data
    measurement_time_index = ones(n_robots,1); % index of last measurement processed
    barcode = 0;
    for i=1:n_robots
        tempIndex=find(Robots{i}.M(:,1)>=start_timestep*sample_time,1,'first');
        if ~isempty(tempIndex)
            measurement_time_index(i) = tempIndex;
        else
            measurement_time_index(i) = n_measurements(i)+1;
        end
    end
    clear tempIndex
    
    % animate robots groundtruth and pose estimates
    for k=start_timestep:end_timestep
        t = k*sample_time;

        if(mod(k,timesteps_per_frame)==0)
            delete(findobj('Type','line')); 
        end

        for i= 1:n_robots        

            x_g(i) = Robots{i}.G(k,2);
            y_g(i) = Robots{i}.G(k,3);
            z_g(i) = Robots{i}.G(k,4);
            
            x_est(i) = Robots{i}.Est(k,2);
            y_est(i) = Robots{i}.Est(k,3);
            z_est(i) = Robots{i}.Est(k,4);

            if(mod(k,timesteps_per_frame)==0)
                x1_g = d_robot*cos(z_g(i)) + x_g(i);
                y1_g = d_robot*sin(z_g(i)) + y_g(i);
                p1_g = x_g(i) - r_robot;
                p2_g = y_g(i) - r_robot;
                set(plotHandles_robot_gt(i),'Position',[p1_g,p2_g,d_robot,d_robot]);
                line([x_g(i) x1_g],[y_g(i) y1_g],'Color','k');
                
                x1_est = d_robot*cos(z_est(i)) + x_est(i);
                y1_est = d_robot*sin(z_est(i)) + y_est(i);
                p1_est = x_est(i) - r_robot;
                p2_est = y_est(i) - r_robot;
                set(plotHandles_robot_est(i),'Position',[p1_est,p2_est,d_robot,d_robot]);
                line([x_est(i) x1_est],[y_est(i) y1_est],'Color','k');
            end

            % plot meaurements of robot i 
            if(draw_measurements)
                while(n_measurements(i) >= measurement_time_index(i) && Robots{i}.M(measurement_time_index(i),1) <= t)
                    measure_id = Robots{i}.M(measurement_time_index(i),2);
                    measure_r = Robots{i}.M(measurement_time_index(i),3);
                    measure_b = Robots{i}.M(measurement_time_index(i),4);
                    landmark_index = find(Barcodes(:,2)==measure_id);
                    if(~isempty(landmark_index))
                        x1 = x(i) + measure_r*cos(measure_b + z(i));
                        y1 = y(i) + measure_r*sin(measure_b + z(i));
                        line([x(i) x1],[y(i) y1],'Color',colour(i,:),'LineWidth',1);
                    else
                        robot_index = find(Barcodes(1:5,2)==measure_id);
                        if(~isempty(robot_index))
                            x1 = x(i) + measure_r*cos(measure_b + z(i));
                            y1 = y(i) + measure_r*sin(measure_b + z(i));
                            line([x(i) x1],[y(i) y1],'Color',colour(i,:),'LineWidth',1);
                        end
                    end
                    measurement_time_index(i) = measurement_time_index(i) + 1; 
                end
            end
        end

        % write time
        if(mod(k,timesteps_per_frame)==0)
            delete(findobj('Type','text'));
            texttime = strcat('k= ',num2str(k,'%5d'), '  t= ',num2str(t,'%5.2f'), '[s]');
            text(1.5,6.5,texttime);
            pause(pause_time_between_frames);
        else
            if(draw_measurements)
                pause(0.001);
            end
        end
    end

    clear Robot x x1 y y1 z r_robot r_landmark p1 p2 max_time measure_id d_robot d_landmark barcode k colour
    clear texttime t measurement_time_index masure_id measure_r measure_b landmark_index robot_index i n_measurements plotHandles* angle

end
