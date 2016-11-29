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
% residual6pntCostFunction is used to find the error with the current
% Kr,Kd,R,andt
%
%Input - Kr - RGB camera calibration matrix
%        Kd - Depth camea calibration matrix
%        R - rotation - R_R_D
%        t - translation - R_t_D
%        projectedCenter_r - RGB camera pixel center of the sphere
%        U_depth - set of points on a sphere
%
%Output - err - residual of the cost function
%
%%

function err = f_residual6pntCostFunction(Kr,Kd,R,t,projectedCenter_r, U_depth)

%Converts pixel center of the sphere to 3D points
X_depth.points = f_depth2XYZ(Kd,U_depth.points);

%Sphere fitting for each sphere
M = f_sphereLinLS(X_depth.points(1:3,:));
centerSphere_hat.center = M(1:3);
radiusSphere_hat = M(4);

% Cost function that calculates the residual
R_P_D = Kr*[R t];
U_temp = R_P_D * [centerSphere_hat.center;1];
U_temp = U_temp/U_temp(3);
err = [projectedCenter_r.points;1] - U_temp;

end

