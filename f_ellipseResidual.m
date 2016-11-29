%%
%Function is used to determine if the points are within a certain threshold
%on the estimated model.
% Uses the equation: r = (a(1-e^2))/(1+e*cos(theta))
% where e = c/a
% where c = sqrt(a^2-b^2)
%
% Input - M - which contains Center - (x,y)
%                           Rx, Ry - (a,b) major, minor axis lengths 
%                           theta  - tilt of the ellipse
%                           phi - angle of each point on the ellipse
%         x - 2D points
%         t - threshold on which points will be considered an outlier or
%         inlier
%
% Output -  inliers -  points on the ellipse within the threshold
%           M
%           outliers - points not on the ellipse
%%


function [inliers, M,outliers,residual] = f_ellipseResidual(M,x,t)


inliers = [];
outliers = [];

z1 = M(1); % x of the center
z2 = M(2); % y of the center
a = M(3); % major/minor axis length
b = M(4); % minor/major axis length
alpha = M(5); % the tilt of the ellipse in radians  

%Rotates the Ellipse back to 0 degrees
R = [cos(alpha) -sin(alpha);
     sin(alpha) cos(alpha)]';
     
%Moves the ellipse back to the origin and at 0 degree
%This is to due to the contraint based on using the atan 2
x_aligned = R*(x - [z1;z2]*ones(1,size(x,2)));

%Computes the Foci of the modeled Ellipse
[F1,AA,c] = f_param2Foci(0,0,a,b,0);

%radius in polar coordinates
r = sqrt((x_aligned(1,:) - F1(1)*ones(1,size(x,2))).^2 + (x_aligned(2,:) - F1(2)*ones(1,size(x,2))).^2);

%eccentricity
ecc = c/a;
%Computes the angle in polar coordinates
t1 = atan2(F1(2) - x_aligned(2,:), F1(1) - x_aligned(1,:));

residual = r - ((a*(1-ecc^2))./(1+ecc*cos(t1)));

indices = abs(residual) <= t^2;
inliers = x(:,indices);
outliers = x(:,logical(1-indices));
end

