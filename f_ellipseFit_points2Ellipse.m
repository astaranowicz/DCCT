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
%Fits an ellipse to a set of data points using Ransac, Linear Least Squares
%and Non-Linear Least Squares
%Input - Points - 2D points
%        threshold - threshold on ransac to pick the inliers
%
%Output - x0,y0 - center of the ellipse
%         a,b  - major, minor axis of the ellipse
%         alpha - tilt of the ellipse
%         phi - angle of each point from the major axis
%         residual - from either the Non-Linear Least Squares or Linear
%         outliers - the points not picked by Ransac
%%

function [x0,y0,a,b,alpha,phi, residual,outliers] = f_ellipseFit_points2Ellipse(Points, threshold)

% toleranceForRotation = 1e-5;
sampleSize = 5; %number of points to describe an ellipse
maxDataTrials = 400; %number of trials of trying to find a non-degenerate model before the code exits
maxTrials = 700; %number of trials before the code exits
inliers = [];
while length(inliers) < sampleSize
    %To pick the best fit to the number of points having no outliners
    [AA, inliers,outliers] = f_ransac_Elp_Sph(Points, @f_ellipseFit2Conic2param, @f_ellipseResidual, sampleSize, threshold,maxDataTrials, maxTrials);
end
%Linear Least Squares to ensure the model is estimated correctly
% [M,par] = f_ellipseLinLS(inliers);
M = f_ellipseFit2Conic2param(inliers);
z1_LS = M(1);
z2_LS = M(2);
a_LS = M(3);
b_LS = M(4);
alpha_LS = M(5);


if a_LS ~= b_LS
    %Generates the points on the ellipse
    Q = [ cos(alpha_LS) -sin(alpha_LS);
        sin(alpha_LS)  cos(alpha_LS)];
    m = size(inliers,2);
    z0 = [z1_LS;z2_LS];
    phi2 = angle([1 1i] * Q' * (inliers - repmat(z0, 1, m)))';
    
    %Non-Linear Least Squares to ensure accuracy in the model
    [NLS,RESIDUAL] = f_ellipseNonLinLS(z1_LS,z2_LS,a_LS,b_LS,alpha_LS,phi2',inliers);
    %Non-LS method minizies all distances to the points and the general form of
    %the ellipse.  The angles to each point is 1:numOfInlierPoints
    sizeOfInliers = length(inliers);
    
    phi = NLS(1:sizeOfInliers);
    alpha = NLS(sizeOfInliers+1);
    a = NLS(sizeOfInliers+2);
    b = NLS(sizeOfInliers+3);
    x0 = NLS(sizeOfInliers+4);
    y0 = NLS(sizeOfInliers+5);
    
    residual = RESIDUAL;
else
    
    alpha = alpha_LS;
    a = a_LS;
    b = b_LS;
    x0 = z1_LS;
    y0 =z2_LS;
    
    residual = 0;
    phi = 0;
    
end

end

