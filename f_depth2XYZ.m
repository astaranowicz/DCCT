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
%Converts the pixel points to the 3D points
%
% Input - Kd - depth camera calibration matrix
%         ud - pixel points with form: [u,v,Z]
%
% Output - XYZ - 3D point in the depth camera
%%

function XYZ = f_depth2XYZ(Kd,ud)

Pd = [inv(Kd) zeros(3,1);
      zeros(1,3)  1];
u = ud(1,:);
v = ud(2,:);
Z = ud(3,:);
%As in eq (6) from the paper
%form: [uZ,vZ,Z];
ud_bar = [u.*Z;v.*Z;Z;ones(1,length(ud(1,:)))];

XYZ = Pd *ud_bar;

end

