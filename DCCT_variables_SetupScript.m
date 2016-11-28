%%
% DEFAULT VALUES  Do not change, unless you do not want the data from the
% paper to run.
%%

DCCT_variables.setOfSpheres = [1:57]; %Number of spheres considered
DCCT_variables.Kd = []; %Depth sensor calibration matrix
DCCT_variables.Kr = [525     0  319.5;
                       0   525  239.5;
                       0     0      1]; %RGB sensor calibration matrix
DCCT_variables.R = []; %Extrinsic calibration parameter, rotation from D to R
DCCT_variables.t = []; %Extrinisc calibration parameter, translation from D to R
DCCT_variables.threshold_3D_Depth = .02; %threshold on DepthMap distance for Sphere fitting in meters
DCCT_variables.threshold_RGB = (.75)*3; %threshold on RGB image for Ellipse fitting in pixels
DCCT_variables.dataPath = 'Data'; %directory where the Data is located
DCCT_variables.RGBImageName = 'rgb'; %name of the RGB Image
DCCT_variables.RGBFileExt = '.png';
DCCT_variables.DepthImageName = 'depth'; %name of the Depth Image
DCCT_variables.depthFileExt = '.png';
DCCT_variables.cannyThreshold_croppedImage = 0.4;
DCCT_variables.cannySigma_croppedImage = 0.15;
DCCT_variables.cannyThreshold_Image = 0.15;
DCCT_variables.cannySigma_Image = 0.25;
DCCT_variables.ProjectSphereCentertolerance = .05; %Tolerance to check if alpha or beta is 0 in projection of sphere center equation
DCCT_variables.threshold_Res6pnt = 10; %pixels - for use in RANSAC Least Squares 6pnt 
DCCT_variables.minDistanceFromCamera = 0.10;%in millimeters for sphere selection in the Depth Map
DCCT_variables.maxDistanceFromCamera = 400; %in millimeters for sphere selection in the Depth Map