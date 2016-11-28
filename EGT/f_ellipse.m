%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox  (EGT)  %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  f_ellipse       Parameterization of an ellipse.
%
%  [x,y] = f_ellipse(A,passo)
%  A = 3by3 symmetric matrix of the ellipse in the essential position.
%  passo = step of parameter
%
%  Descr: 
%  -----  This function returns the parameterization of the ellipse in the
%         essential position associated with the matrix. The known form of
%         an ellipse equation in the essential position is 
%         x^2/a^2+y^2/b^2=1
%         and the parameterization returned by the function is
%         [x,y] = [a*cos(t),-b*sin(t)] where 0<=t<=(2*pi)
%

function [u,v] = f_ellipse(A,passo,ampiezza);
if nargin==1,
     passo=pi/360;
 elseif nargin==2,
     ampiezza=2*pi;
 end;
 ampiezza=2*pi;
 passo=pi/360;
a = sqrt(-A(3,3)/A(1,1));
b = sqrt(-A(3,3)/A(2,2));
t = 0:passo:ampiezza;
u = a*cos(t);
v = -b*sin(t);


