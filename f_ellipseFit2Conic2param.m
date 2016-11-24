%%
%ellipseFit2Conic2param is a linear least squares estimation of an
%ellipse,which moves from the vector of parameters to the parametric form of
%the ellipse
%Input -  Points - 2D points that are used to fit an ellipse
%       
%Output - param - vector of parametric parameters containing:
%                 [x0,y0,a,b,alpha] ~ [center, semiaxis a, semiaxis b, tilt
%                 of the ellipse]
%%

function param = f_ellipseFit2Conic2param(Points)

% A = [a b c d e f]' is the vector of parameters of the fitting ellipse:
% ax^2 + bxy + cy^2 +dx + ey + f = 0
%Ensures that the points are stacked column-wise, not row-wise
if length(Points) > 2
   Points = Points'; 
end
% Linear LeastSquares ellipse fit
%Based on "Numerically Stable Direct Least Squares Fitting Of Ellipses" by Halir and Flusser
%By Fitzgibbon
%  A= f_ellipseLinLS(Points');
%Based on "Numerically Stable Direct Least Squares Fitting Of Ellipses" by Halir and Flusser
%By another person
A = f_EllipseDirectFit(Points);
check = find(A == 0);
if length(check) == 6
    param = [0,0,0,0,0];
    return
end
%Normalization to force the last element to be +/-1 
% A = A/abs(A(end));
% A = A/A(1);
%A = f_roundn(A,-16);

if isreal(A) == 0
    param = [0,0,0,0,0];
    return
end

%Input - the vector of parameters from the LS 
%      - toleranceForImConic for if-statements
%Output - parametric form parameters
% toleranceForRotation = 1e-5;
[AAA,a,b,R,t] = imconic([A(1) A(2) A(3) A(4) A(5) A(6)]);
%R is the Rotation matrix which contains the tilt of the ellipse
%   alpha is the angle associated with the Rotation matrix
%   atan2 is used to ensure the angle is between 0-180
%   atan is not used due to the ratio,  It can not discriminate between
%   the angles of 45/-45
alpha = atan2(R(2,1),R(2,2));


%Vector of parameters containing the center, semiaxis a, semiaxis b, and
%tilt
param = [t(1),t(2),a,b,alpha];

end
