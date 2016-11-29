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
%Converts the Paramertic form of the Ellipse to the Conic
%
%Taken from: Charles F. Van Loan, Dept. C.S, Cornell University
%            "Using the Ellipse to Fit and Enclose Data Points"
%            http://www.cs.cornell.edu/cv/OtherPdf/Ellipse.pdf
%
%Input- [x0,y0] - center of ellipse
%       a,b  - major, minor axis of the ellipse
%       alpha - tilt of the ellipse
%
%Output - Conic - the conic form of the ellipse Ax^2+Bxy+Cy^2+Dx+Ey+F=0
%
%
%Note: special case is when alpha is pi/4,  means that s = c
%%
function Conic = f_param2Conic_Ellipse(x0,y0,a,b,alpha)

s = sin(alpha);
c = cos(alpha);

A = (b*c)^2 + (a*s)^2;
B = -2*c*s*(a^2-b^2);
C = (b*s)^2 + (a*c)^2;
D = -2*A*x0 - y0*B;
E = -2*C*y0 - x0*B;
F = -(a*b)^2 + A*x0^2+ B*x0*y0 + C *y0^2;

% Normalization to ensure the last element is 1 and positive
v=[A B C D E F]';
Conic = [ v(1)   v(2)/2  v(4)/2;
          v(2)/2 v(3)    v(5)/2;
          v(4)/2 v(5)/2  v(6)];

end

