%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox  (EGT)  %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  f_hyperbola       Parameterization of a hyperbola.
%
%  [x,y] = f_hyperbola(A,passo)
%  A = 3by3 symmetric matrix of the hyperbola in the essential position.
%  passo = step of parameter
%
%  Descr: 
%  -----  This function returns the parameterization of the positive
%         sheet of the hyperbola in the essential position, associated 
%         with the matrix, to plot the hyperbola it has to be considered
%         also the negative sheet. 
%         The essential equation of a hyperbola is
%         x^2/a^2-y^2/b^2=1
%         and the parameterization returned by the function is
%         [x,y] = [t,(+-)a*sqrt(1+t^2/b^2)] 
%         where the parameter t is bounded by the image resolution.
%
%  Example:
%       A = [1 0 0; 0 -1 0; 0 0 -1];
%       [x,y] = f_hyperbola(A,passo);
%       figure(1),hold on, 
%       plot(x,y,'b');
%       plot(x,-y,'b');
%       axis equal;
%

function [u,v] = f_hyperbola(A,passo,ampiezza);
if nargin==1,
     passo=1/100;
 elseif nargin==2,
     ampiezza=10;
 end;
 
t = -ampiezza:passo:ampiezza;%era 10
a = sqrt(-A(3,3)/A(1,1));
b = sqrt(A(3,3)/A(2,2));
v = t;
u = a*sqrt(1+t.^2/b^2);

% if A(1,1)*A(3,3)>0
%    a = sqrt(-A(3,3)/A(2,2));
%    b = sqrt(A(3,3)/A(1,1));
%    u = t; 
%    v = a*sqrt(1+t.^2/b^2);
% %%%%%%%% vneg = - a*sqrt(1+t.^2/b^2);
% else
%    a = sqrt(-A(3,3)/A(1,1));
%    b = sqrt(A(3,3)/A(2,2));
%    u = -50:.01:50;
%    v = b*sqrt(-1+t.^2/a^2);
% end

