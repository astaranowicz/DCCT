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
% DepthCal_points uses the linear least squares method to estimate the
% depth camera calibration.  The DLT (6point) algorithm used is a modified
% version that uses  || |Kr*u_r|_x - [R t]K_depth*u_bar_d|| = 0
% where u_r are the center of the sphere in pixel from the rgb camera
% and u_bar_d is the center of the sphere in pixel from the depth map
%
% Input - pixelCenter - center of sphere from RGB camera in pixel
%         centersOfSphere - center of sphere from Depth Map in pixel
%         Kr - camera calibration of the RGB camera
%         threshold_Res6pnt - threshold for the residuals used to pick the
%                             inliers
%
% Output - Kd_est - estimated depth camera calibration
%          R_est - estimated rotation from R_R_D
%          t_est - estimated translation from R_t_D
%          Res - residual calulated for each sphere based on the estimated
%                 model
%          inlier_indicies - indicies of the spheres that where used to
%                             find the best estimated model from RANSAC
%
%%

function [Kd_est, R_est, t_est,Res,inlier_indicies] = f_DepthCal_points(pixelCenter,centersOfSphere,Kr,threshold_Res6pnt)

sampleSize = 6; %minimal number of points needed for the DLT algorithm (6point)
maxDataTrials = 500; %number of trials of trying to find a non-degenerate model before the code exits
maxTrials = 1000; %number of trials before the code exits

%RANSAC to find the best spheres to estimate the model
[AA, inlier_indicies, outlier_indicies] = f_ransac_6Pnt(pixelCenter, centersOfSphere,Kr, ...
    @f_6pntalgorithmKd, @f_residual6Point, sampleSize, threshold_Res6pnt,maxDataTrials, maxTrials);

%Linear Least Squares (DLT) to find the estimated Kd from the set of inliers
M= f_6pntalgorithmKd(pixelCenter(:,inlier_indicies), centersOfSphere(:,inlier_indicies), Kr);
%Calculate the Residual based on the inliers and best-fit model
[AA, BB, CC, Res] = f_residual6Point(M,pixelCenter(:,inlier_indicies), centersOfSphere(:,inlier_indicies),Kr, threshold_Res6pnt);

Kd_est =M.Kd;
R_est = M.R;
t_est = M.t;

end

