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
%Compute the Euler Angles
%Input - R - rotation matrix
%
%Output -
%          row - z axis
%          pitch - y axis
%          yaw - x axis
%%


function [roll,pitch,yaw] = f_rotationMat2EulerAngles(R)

%From Siciliano Book
%The known parameters for the Rotation Matrix
%pitch range = (-pi/2, pi/2)

roll_est = atan2(R(2,1), R(1,1));
pitch_est = atan2(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));
yaw_est = atan2(R(3,2), R(3,3));
%pitch range = (pi/2, 3pi/2)
roll2_est = atan2(-R(2,1), -R(1,1));
pitch2_est = atan2(-R(3,1), -sqrt(R(3,2)^2 + R(3,3)^2));
yaw2_est = atan2(-R(3,2), -R(3,3));

if(pitch_est > -pi/2 && pitch_est < pi/2),
    roll = roll_est;
    pitch = pitch_est;
    yaw = yaw_est;
else
    roll = roll2_est;
    pitch = pitch2_est;
    yaw = yaw2_est;
end


end

