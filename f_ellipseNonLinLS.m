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

