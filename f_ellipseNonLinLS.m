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
%Fits an ellipse using non-linear least squares
%Input -  Initial z, a, b, alpha, and phi from linear least squares
%
%Output - Estimated z, a, b, alpha, and phi from non-linear least squares
%%

function [NLS,RESIDUAL] = f_ellipseNonLinLS(z1,z2,a,b,alpha,phi,Points)

global PointsForNLSEllipse

%Points that lie on the ellipse
PointsForNLSEllipse = Points;
%Variables that will be minimized- Initial estimates
%Order matters with this set up: From Gander_BIT94
Z = [phi'; alpha; a; b ;z1; z2]';

options = optimset('Algorithm','levenberg-marquardt','TolX',1e-5,'TolFun',1e-5,  ... 
                      'Jacobian','on','MaxIter',600,'MaxFunEvals',600,'Display','off');

lb = [];
ub = [];

[NLS,~,RESIDUAL,~,~,~,JACOBIAN] = lsqnonlin(@f_minEllipse,Z,lb,ub,options);


end

