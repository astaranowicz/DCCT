%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C) 2016 Aaron Staranowicz and Gian Luca Mariottini
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%Residual for the DLT (6point) algorithm cost function LS
%
% Input - M - Struct containing: Kd, R, t
%                   Kd - depth camera calibration matrix
%                   R - R_R_D - rotation from RGB camera to Depth Camera
%                   t - R_t_D - translation from RGB camera to Depth Camera
%       pixelCenter - center of sphere in pixel from RGB camera
%       centersOfSphere - center of sphere in pixel from Depth Map
%       Kr - RGB camera calibration matrix
%       threshold - threshold to select the inliers
%
% Output - inliers - inliers that met the threshold
%          M - struct containing: Kd,R,t
%          outliers - outliers that did not met the threshold
%          residual - residual based on the calculations
%
%%

function [inliers, M, outliers, residual] = f_residual6Point(M,pixelCenter,centersOfSphere,Kr, threshold)

% R_t_D
t = M.t;
% R_R_D
R = M.R;
% depth camera calibration
Kd = M.Kd;
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
X_depth = f_depth2XYZ(Kd,centersOfSphere);
%Uses the original 6pnt algorithm: ||u_r - Kr[R t]X_d||
for i = 1:length(centersOfSphere)
    X_depth_temp = temp_KrRt * X_depth(:,i);
    U_depth(:,i) = X_depth_temp / X_depth_temp(3);
    residual(i) = norm([pixelCenter(:,i);1] - U_depth(:,i));
end

indices =  residual <= threshold;
inliers = find(indices);
outliers = find(logical(1-indices));

end

