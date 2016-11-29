%%
%Residual for the DLT (6point) algorithm cost function LS
%
% Input - M - Struct containing: Kd, Dd, R, t
%                   Kd - depth camera calibration matrix
%                   Dd - depth camera distortion vector
%                   R - R_R_D - rotation from RGB camera to Depth Camera
%                   t - R_t_D - translation from RGB camera to Depth Camera
%       pixelCenter - center of sphere in pixel from RGB camera
%       centersOfSphere - center of sphere in pixel from Depth Map
%       Kr - RGB camera calibration matrix
%       Dr - RGB camera distorion vector
%       threshold - threshold to select the inliers
%
% Output - inliers - inliers that met the threshold
%          M - struct containing: Kd,Dd,R,t
%          outliers - outliers that did not met the threshold
%          residual - residual based on the calculations
%
%%

function [inliers, M, outliers, residual] = f_residual6Point(M, pixelCenter, centersOfSphere, Kr, Dr, threshold)

% R_t_D
t = M.t;
% R_R_D
R = M.R;
% depth camera calibration
Kd = M.Kd;
Dd = M.Dd;
%Checks if the estimated model that was passed was the special case of
%being ill-coniditioned from the function: f_6pntalgorithmKd
if Kd(1,1) == 0
    inliers = [];
    outliers = 1:length(centersOfSphere);
    residual = [];
    return;
end

temp_KrRt = Kr*[R t];
%Converts the pixels from the Depth Map to 3D points
X_depth = f_depth2XYZ(Kd,Dd,centersOfSphere);
%Uses the original 6pnt algorithm: ||u_r - Kr[R t]X_d||
for i = 1:length(centersOfSphere)
    XYZ = X_depth(1:3,i);
    XYZ = f_undistort(XYZ, Dr);
    X_depth_temp = temp_KrRt * [XYZ;1];
    U_depth(:,i) = X_depth_temp / X_depth_temp(3);
    residual(i) = norm([pixelCenter(:,i);1] - U_depth(:,i));
end

indices =  residual <= threshold;
inliers = find(indices);
outliers = find(logical(1-indices));

end

