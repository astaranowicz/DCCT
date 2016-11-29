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
%Function is used to determine if the points are within a certain threshold
%on the estimated model.
%Residual for a sphere follows this equation:
% d_i = ||Xi - X0|| - r
%
%Input - M - (x,y,z) center and (r) radius
%        x - points of the sphere
%        t - threshold to check if the points lie on the sphere or not
%
%Output - inliers - the points on the sphere with a certain threshold
%         M - (x,y,z) center and (r) radius
%         outliers - the points not on the sphere
%%

function [inliers, M, outliers,indices] = f_sphereResidual(M,x,t)

inliers = [];
outliers = [];

x0 = M(1); % x center
y0 = M(2); % y center
z0 = M(3); % z center
radius_i = M(4);  % radius of sphere

x2 = x- [x0; y0; z0] * ones(1, size(x,2));

dist = f_roundn(sqrt(sum(x2 .* x2)) - radius_i,-10);

indices = dist.^2 <= t^2;
inliers = x(:,indices);
outliers = x(:,logical(1-indices));

end




