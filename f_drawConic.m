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
% Draws a conic using imconic's drawing function
%
% Input - Conic - the conic that will be drawn
%         fighandle - handle of the figure to draw the ellipse
%         color -  color of the ellipse edge
%%

function f_drawConic(Conic,fighandle, color)
A = Conic(1,1);
B = Conic(1,2)*2;
C = Conic(2,2);
D = Conic(1,3)*2;
E = Conic(2,3)*2;
F = Conic(3,3);

imconic([A B C D E F],fighandle,color,[]);

end


