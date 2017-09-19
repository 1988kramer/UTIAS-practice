addpath ../common;
%addpath ../MRCLAM_Dataset1;

n_robots = 1;
[Barcodes, Landmark_Groundtruth, Robots] = loadMRCLAMdataSet(n_robots);

% rest of the code goes here