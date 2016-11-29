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
%Calculates the Projection of the Sphere in the image from the estimated
%center of the ellipse
%
%Input - x0,y0 - center of the ellipse
%        a,b  - major, minor axis of the ellipse
%        alpha - tilt of the ellipse
%        K_camera - the calibration matrix for the camera
%                       form: [ fku 0  u0;
%                                0 fkv v0;
%                                0  0   1];
%
%Output - sphere_center - (u,v) of the projected sphere center from the
%                               ellipse center
%%

function sphere_center = f_projectionSphere(x0, y0, a, b, alpha, K_camera,tolerance)

%To find the ratio of R/C to be used to find the sphere projected center
gamma_r = f_findRatioR_C(x0, y0, a, b, alpha, K_camera,tolerance);

%The center of the projected sphere in the camera frame
u = x0 * (1 - gamma_r(1)) + K_camera(1,3) * gamma_r(1);
v = y0 * (1 - gamma_r(2)) + K_camera(2,3) * gamma_r(2);

sphere_center = [u;v];
%  sphere_center(1) =  (1 - r_c) * (x0 - K_camera(1,3));
%  sphere_center(2) =  (1 - r_c) * (y0 - K_camera(2,3));

end
