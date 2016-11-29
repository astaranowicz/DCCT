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
%Changes from the parameteric form to the foci of the ellipse
%
%
%Input - z1,z2 - center of the ellipse
%        a,b - major/minor axis length
%        alpha - tilt of the ellipse
%
%Output - F1 - focus of the first 
%         F2 - 2nd focus
%         c - distance between a and b
%%

function [F1,F2,c] = f_param2Foci(z1,z2,a,b,alpha)



c = sqrt(abs(a^2 - b^2));

F1 = [z1 - (cos(alpha)*c);
      z2 - (sin(alpha)*c)];

F2 = [z1 + (cos(alpha)*c);
      z2 + (sin(alpha)*c)];

end

