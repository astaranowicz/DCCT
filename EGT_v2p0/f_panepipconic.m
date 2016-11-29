%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%
%  %  Epipolar Geometry Toolbox v1.3 (EGT)  %%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  f_panepipconic(nvers,K,color,a,b,quadric,passo,ampiezza);
%
% Description:
% -----------
%  Plots the epipolar conics in the CCD image plane (pixels)
%
% Syntax:
% ------
%  nvers = normal versor to the epipolar plane. Note that, being t the
%          translation centered at the desired and pointing toward the current
%          view, then nvers=E'*Xhd   and nvers'=E*Xhc where E=f_skew(t)*R
%          and Xhc(d) is the mirror projection in the current (desired) camera 
%          of the scene point P.
%  K = calibration matrix with the inner parameters of the CCD camera.
%  color = character strings which defines line type and colour used to plot the 
%          epipolar conic.
%  a,b = Mirror parameters.
%  quadric = 1 or 2 for hyperboloid or paraboloid, respectively.
%  passo, ampiezza= parameters for plot of conics parameterizations.
% 
% Authors:
%    Gian Luca Mariottini
% Thanks to:
%    Nicola Pisu
%
% Last update:
%    December, 2005
%
function [A22,flag] = f_panepipconic(nvers,K,color,a,b,quadric,passo,ampiezza);

if nargin==1
    K = eye(3);
    color = 'b';
    a = 3;
    b = 1;
    passo=1/100;
    ampiezza=10;
elseif nargin==2
    color = 'b';
    a = 3;
    b = 1;
    passo=1/100;
    ampiezza=10;
elseif nargin==3
    color = 'b';
    a = 3;
    b = 1;
    passo=1/100;
    ampiezza=10;
elseif nargin==4
    a = 3;
    b = 1;
    passo=1/100;
    ampiezza=10;
elseif nargin==5
    b = 1;
    passo=1/100;
    ampiezza=10;
    norma=1;
elseif nargin==6
    passo=1/100;
    ampiezza=10;
elseif nargin==7;
    ampiezza=10;
elseif nargin>8,    
    display('  EGT error: too many inputs in f_panepipconic')
end;


% test to verify if the nvers is in homogeneous notation %
if length(nvers)==4
    nvers([1:3]) = nvers([1:3])/(nvers(3));
end
    
% projection of the conic in the CCD image plane %
Rc = eye(3);

p = nvers(1);
q = nvers(2);
s = nvers(3);
if quadric==2,
A = [ s,      0      ,    a*p    ;           
      0,      s      ,    a*q    ;
     a*p,    a*q     ,    -a^2*s];
elseif quadric==1
e=sqrt(a^2+b^2);
A= [ -4*s^2*a^2*e^2+p^2*b^4  ,      p*q*b^4          ,    p*s*b^2*(-2*e^2+b^2) ;           
       p*q*b^4               , -4*s^2*a^2*e^2+q^2*b^4,    q*s*b^2*(-2*e^2+b^2) ;
       p*s*b^2*(b^2-2*e^2)   , q*s*b^2*(b^2-2*e^2)   ,              s^2*b^4   ];
end;   
A2 = inv(K)'*Rc*A*Rc'*inv(K);
A22=A2/max(max(A2)); %Altrimenti sbaglia un po'...non so bene perche'!

flag = f_conics(A22,color,passo,ampiezza);
