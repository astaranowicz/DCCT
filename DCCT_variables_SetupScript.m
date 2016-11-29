%%
% DEFAULT VALUES  Do not change, unless you do not want the data from the
% paper to run.
%%

DCCT_variables.setOfSpheres = [1:25]; %Number of spheres considered
DCCT_variables.Kr = [ 506.522971745157750  0  329.706020747618250;
                        0   506.786072854127440  265.893095027934810;
                        0   0   1]; %RGB sensor calibration matrix
DCCT_variables.Kd = []; %Depth sensor calibration matrix
DCCT_variables.R = []; %Extrinsic calibration parameter, rotation from D to R
DCCT_variables.t = []; %Extrinisc calibration parameter, translation from D to R
DCCT_variables.threshold_3D_Depth = .02; %threshold on DepthMap distance for Sphere fitting in meters
DCCT_variables.threshold_RGB = (.75)*3; %threshold on RGB image for Ellipse fitting in pixels
DCCT_variables.dirExt = 'Data/'; %directory where the Data is located
DCCT_variables.RGBImageName = 'RGBImage_ball_'; %name of the RGB Image
DCCT_variables.DepthImageName = 'DepthImage_ball_'; %name of the Depth Image
DCCT_variables.fileExt = '.jpg';
DCCT_variables.DepthMapName = 'ball_'; %name of the Depth map
DCCT_variables.fileExtDM = '.txt';
DCCT_variables.ProjectSphereCentertolerance = .05; %Tolerance to check if alpha or beta is 0 in projection of sphere center equation
DCCT_variables.threshold_Res6pnt = 10; %pixels - for use in RANSAC Least Squares 6pnt 
DCCT_variables.minDistanceFromCamera = 0.10;%in millimeters for sphere selection in the Depth Map
DCCT_variables.maxDistanceFromCamera = 400; %in millimeters for sphere selection in the Depth Map