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
%Fits a sphere to a set of data points using Ransac, Linear Least Squares
%and Non-Linear Least Squares
%Input - Points - 3D points
%        threshold - threshold on ransac to select the inliers
%
%Output - center - center of the sphere (x,y,z)
%         radius - radius of the sphere
%         residual - residual from Non-Linear Least Squares
%         inliers - points on the sphere within the threshold
%         outliers - points not on the sphere
%
%%

function [center, radius, residual, inliers, outliers, indicesReal] = f_sphereFit_points2Sphere(Points,threshold)

sampleSize = 4; %number of points to describe an ellipse
maxDataTrials = 600; %number of trials of trying to find a non-degenerate model before the code exits
maxTrials = 1000; %number of trials before the code exits

%To pick the best fit to the number of points having no outliners
[AA, inliers,outliers,logicalIndices] = f_ransac_Elp_Sph(Points, @f_sphereLinLS, @f_sphereResidual, sampleSize, threshold,maxDataTrials, maxTrials);

indicesReal = find(logicalIndices);
%Linear Least Squares to find the best-fit model with inliers
M = f_sphereLinLS(inliers);

center_LS = M(1:3);
radius_LS = M(4);

%Non-Linear Least Squares to minimize the residual
%[center_NLS,radius_NLS,RESIDUAL_NLS] = f_sphereNonLinLS(center_LS, radius_LS,inliers);

center = center_LS;
radius = radius_LS;
residual = 0;
end
