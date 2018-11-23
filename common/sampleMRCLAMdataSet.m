% UTIAS Multi-Robot Cooperative Localization and Mapping Dataset
% produced by Keith Leung (keith.leung@robotics.utias.utoronto.ca) 2009
% Matlab script animateMRCLAMdataSet.m
% Description: This scripts samples the dataset at fixed intervals 
% (default is 0.02s). Odometry data is interpolated using the recorded time. % Measurements are rounded to the nearest timestep. 
% Run this script after loadMRCLAMdataSet.m
function [Robots, timesteps] = sampleMRCLAMdataSet(Robots, sample_time)

    % sample_time = 0.02;
    
    n_robots = size(Robots, 1);
    
    min_time = Robots{1}.G(1,1);
    max_time = Robots{1}.G(end,1);
    for n=2:n_robots
       min_time = min(min_time, Robots{n}.G(1,1));
       max_time = max(max_time, Robots{n}.G(end,1));
    end
    for n=1:n_robots
        Robots{n}.G(:,1) = Robots{n}.G(:,1) - min_time;
        Robots{n}.M(:,1) = Robots{n}.M(:,1) - min_time;
        Robots{n}.O(:,1) = Robots{n}.O(:,1) - min_time;
    end
    max_time = max_time - min_time;
    timesteps = floor(max_time/sample_time)+1;

    oldData = 0;
    for n = 1:n_robots
        oldData = Robots{n}.G;

        k = 0;
        t = 0;
        i = 1;
        p = 0;

        [nr,nc] = size(oldData);
        newData = zeros(timesteps,nc);
        while(t <= max_time)
            newData(k+1,1) = t;     
            while(oldData(i,1) <= t);        
                if(i==nr)
                    break;
                end
                i = i + 1;
            end
            if(i == 1 || i == nr)
                newData(k+1,2:end) = 0;
            else
                p = (t - oldData(i-1,1))/(oldData(i,1) - oldData(i-1,1));
                if(nc == 8) % i.e. ground truth data
                    sc = 3;
                    newData(k+1,2) = oldData(i,2); % keep id number 
                else
                    sc = 2;
                end
                for c = sc:nc
                    if(nc==8 && c>=6)
                        d = oldData(i,c) - oldData(i-1,c);
                        if d > pi
                            d = d - 2*pi;
                        elseif d < -pi
                            d = d + 2*pi;
                        end
                        newData(k+1,c) = p*d + oldData(i-1,c);
                    else
                        newData(k+1,c) = p*(oldData(i,c) - oldData(i-1,c)) + oldData(i-1,c);
                    end
                end
            end
            k = k + 1;
            t = t + sample_time;
        end

        Robots{n}.G = newData ;
    end
    
    oldData = 0;
    for n = 1:n_robots
        oldData = Robots{n}.O;

        k = 0;
        t = 0;
        i = 1;
        p = 0;

        [nr,nc] = size(oldData);
        newData = zeros(timesteps,nc);
        while(t <= max_time)
            newData(k+1,1) = t;     
            while(oldData(i,1) <= t);        
                if(i==nr)
                    break;
                end
                i = i + 1;
            end
            if(i == 1 || i == nr)
                newData(k+1,2:end) = oldData(i,2:end);
            else
                p = (t - oldData(i-1,1))/(oldData(i,1) - oldData(i-1,1));
                if(nc == 8) % i.e. ground truth data
                    sc = 3;
                    newData(k+1,2) = oldData(i,2); % keep id number 
                else
                    sc = 2;
                end
                for c = sc:nc
                    if(nc==8 && c>=6)
                        d = oldData(i,c) - oldData(i-1,c);
                        if d > pi
                            d = d - 2*pi;
                        elseif d < -pi
                            d = d + 2*pi;
                        end
                        newData(k+1,c) = p*d + oldData(i-1,c);
                    else
                        newData(k+1,c) = p*(oldData(i,c) - oldData(i-1,c)) + oldData(i-1,c);
                    end
                end
            end
            k = k + 1;
            t = t + sample_time;
        end

        Robots{n}.O = newData ;
    end
    
    oldData = 0;
    for n = 1:n_robots
        oldData = Robots{n}.M;
        newData=oldData;
        for i = 1:length(oldData)
            newData(i,1) = floor(oldData(i,1)/sample_time + 0.5)*sample_time; 
        end
        Robots{n}.M = newData ;
    end

    clear min_time oldData newData nr nc n p sc t k c i d array_names name;
end

