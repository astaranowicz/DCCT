%%
%Fits an ellipse using the algebraic equation
%Used from Numerically Stable Direct Least Squares Fitting Of Ellipses by
%Halir and Flusser
%
%Inputs - Points which contains the (X,Y) of the data set
%
%Output - M  which contains Center - (x,y)
%                           Rx, Ry - (a,b) major, minor axis lengths 
%                           theta  - tilt of the ellipse
%                           phi - angle of each point on the ellipse
%%

function par = f_ellipseLinLS(Points)

X = Points(1,:);
Y = Points(2,:);

% normalize data
mx = mean(X);
my = mean(Y);
sx = (max(X)-min(X))/2;
sy = (max(Y)-min(Y))/2; 

x = (X-mx)/sx;
y = (Y-my)/sy;

% Force to column vectors
x = x(:);
y = y(:);

% Build design matrix
D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ];

% Build scatter matrix
S = D'*D;

% Build 6x6 constraint matrix
C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;

% Solve eigensystem

% New way, numerically stabler in C [gevec, geval] = eig(S,C);

% Break into blocks
tmpA = S(1:3,1:3); 
tmpB = S(1:3,4:6); 
tmpC = S(4:6,4:6); 
tmpD = C(1:3,1:3);
tmpE = inv(tmpC)*tmpB';
[evec_x, eval_x] = eig(inv(tmpD) * (tmpA - tmpB*tmpE));

% Find the positive (as det(tmpD) < 0) eigenvalue
I = find(real(diag(eval_x)) < 1e-8 & ~isinf(diag(eval_x)));

% Extract eigenvector corresponding to negative eigenvalue
A = real(evec_x(:,I));

% Recover the bottom half...
evec_y = -tmpE * A;
A = [A; evec_y];


% unnormalize
par = [
  A(1)*sy*sy,   ...
      A(2)*sx*sy,   ...
      A(3)*sx*sx,   ...
      -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
      -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
      A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
      - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
      + A(6)*sx*sx*sy*sy   ...
      ]';

% % Convert to geometric radii, and centers
% 
% thetarad = 0.5*atan2(par(2),par(1) - par(3));
% %thetarad = 0.5*atan(par(2)/(par(1) - par(3))); %GLM:uncomment the above line
% display(['theta_before_if = ',num2str(thetarad*180/pi)]);
% %thetarad = -0.5*atan2(par(2),par(3)- par(1));
% cost = cos(thetarad);
% sint = sin(thetarad);
% sin_squared = sint.*sint;
% cos_squared = cost.*cost;
% cos_sin = sint .* cost;
% 
% Ao = par(6);
% Au =   par(4) .* cost + par(5) .* sint;
% Av = - par(4) .* sint + par(5) .* cost;
% Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
% Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;
% 
% % ROTATED = [Ao Au Av Auu Avv]
% 
% tuCentre = - Au./(2.*Auu);
% tvCentre = - Av./(2.*Avv);
% wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;
% 
% uCentre = tuCentre .* cost - tvCentre .* sint;
% vCentre = tuCentre .* sint + tvCentre .* cost;
% 
% Ru = -wCentre./Auu;
% Rv = -wCentre./Avv;
% 
% Ru = sqrt(abs(Ru)).*sign(Ru);
% Rv = sqrt(abs(Rv)).*sign(Rv);
% 
% z1 = uCentre;
% z2 = vCentre;
% a = Ru;
% b = Rv;
% alpha = thetarad;
% 
% %Changes the alpha to be between 0-180 
% % if sign(roundn(alpha,-3)) < 0,
% %     alpha  = alpha+pi/2;
% % 
% %     a_temp = a;
% %     b_temp = b;
% %     
% %     a = b_temp;
% %     b = a_temp;
% % 
% % end
% %Imposes the constraint that A has to be greater than B
% if a < b,
%     alpha  = alpha+(pi/2);%GLM This was "-"
% 
%     a_temp = a;
%     b_temp = b;
%     
%     a = b_temp;
%     b = a_temp;
%     
% end
% 
% M = [z1;z2;a;b;alpha];



end



