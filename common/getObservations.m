function [z, index] = getObservations(Robots, robot_num, t, index, codeDict)
    % build vector of features observed at current time
    z = zeros(3,1);
    while (Robots{robot_num}.M(index, 1) - t < .005) && (index < size(Robots{robot_num}.M,1))
        barcode = Robots{robot_num}.M(index,2);
        landmarkID = 0;
        if (codeDict.isKey(barcode))
            landmarkID = codeDict(barcode);
        else
            disp('key not found');
        end
        if landmarkID > 5 && landmarkID < 21
            range = Robots{robot_num}.M(index, 3);
            bearing = Robots{robot_num}.M(index, 4);
            if uint8(z(3)) == 0
                z = [range;
                     bearing;
                     landmarkID - 5];
            else
                newZ = [range;
                        bearing;
                        landmarkID - 5];
                z = [z newZ];
            end
        end
        index = index + 1;
    end
end