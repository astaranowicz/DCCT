%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox  (EGT)  %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  f_parabola      Parameterization of a parabola.
%
%  [x,y] = f_parabola(A,passo)
%  A = 3by3 symmetric matrix of the parabola in the essential position.
%  passo = step of parameter
%
%  Descr: 
%  -----  This function returns the parameterization of the parabola in the
%         essential position associated with the matrix. 
%

function [u,v] = f_parabola(A,passo,ampiezza);
 if nargin==1,
     passo=1/100;
     ampiezza=10;
 elseif nargin==2,
     ampiezza=10;
 end;
 
t = -ampiezza:passo:ampiezza; % cambiare (limiti dell'immagine)

if A(2,2)==0
    u = t;
    v = -(A(1,1)/2*A(2,3))*t.^2;
elseif A(1,1)==0
    v = t;
    u = -(A(2,2)/2*A(1,3))*t.^2;
end;
    