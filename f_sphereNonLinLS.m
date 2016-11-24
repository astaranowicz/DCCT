%%
%Estimates the x,y,z and radius of a sphere using Non-Linear Least Squares
%
% Inputs - X0- Initial Estimate of the center [x,y,z] and radius r
%
%
% Outputs- Estimated Center [x,y,z] and Radius r
%%

function [center_NLS,radius_NLS,RESIDUAL_i] = f_sphereNonLinLS(center_hat, radius_hat,Depth_points)

global PointsForNLS
% Points on the sphere
PointsForNLS = Depth_points;
% The variables to be minizied
Z = [center_hat(1); center_hat(2); center_hat(3); radius_hat];

options = optimset('Algorithm','levenberg-marquardt','TolX',0.00001,'TolFun',0.00001,  ... 
                      'Jacobian','off','MaxIter',600,'MaxFunEvals',600,'Display','off');

[NLS,RESNORM,RESIDUAL_i] = lsqnonlin(@f_minSphere,Z,[],[],options);


center_NLS = [NLS(1);NLS(2);NLS(3)];
radius_NLS = NLS(4);

end



