%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox  (EGT)  %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [Xhmir] = f_daqaXhmir(Q,K,a,b,quadric)
%
%  Reprojects image points to the central catadioptric camera mirror.
%   
%  DESCR: 
%  ------  
%         This function reprojects the image points Q (given in pixel
%         coordinates and obtained through a central catadioptric camera) 
%         into mirror points Xhmir, and returns the coordinates of these 
%         ones expressed in the mirror coordinate system.
%
% SYNTAX:
% ------
%  Q = 3-by-n matrix (if the points are in homogeneous coordinates) or 
%      2-by-n matrix (no homogeneous coordinates) with image points (pixel 
%      coordinates) to reproject. 
%  K = camera calibration matrix with the internal parameters of the camera.
%      K is a 3by3 matrix. For default K is equal to identity.
%  a,b = mirror parameters (scalar). If no specified, default
%        values are a = 3 cm and b = 1 cm.
%  quadric = 1 or 2 for hyperbolic or parabolic mirror, respectively.
%
%Authors:
%  Gian Luca Mariottini 
%  Eleonora Alunno, 
%
% Last update:
%  December 2005.
%

function [Xhmir] = f_daqaXhmir(Q,K,a,b,quadric);

if nargin==1
    K = eye(3);
    a = 3;
    b = 1;
    quadric=1;
elseif nargin==2
    a = 3;
    b = 1;
    quadric=1;
elseif nargin==3
    b = 1;
    quadric=1;
elseif nargin==4,
    quadric=1;    
elseif nargin>5
    display('  EGT error: too many inputs in f_daqaXhmir');
end

% control to verify if the vectors inside the matrix Q 
% are homogeneous 
if length(Q(:,1))==2
    Q = [Q;ones(1,length(Q(1,:)))];
end

e = sqrt(a^2+b^2);
Rcam = eye(3);  % Rcam = Rcam2mir
tcam = [0 0 -2*sqrt(a^2+b^2)]'; % tcam = tcam2mir
[m,n] = size(Q);

if quadric==1, %Hyperbolic mirror
for i=1:n,
    q = Q(:,i);
    v = Rcam'*inv(K)*q;
    if (b^2*v(3)^2-a^2*v(1)^2-a^2*v(2)^2)==0
        Fh = (b^2)/(2*e*v(3));
    else
        Fh = b^2*(e*v(3)+a*norm(v))/(b^2*v(3)^2-a^2*v(1)^2-a^2*v(2)^2);
    end
    Xhmir(:,i) = Fh*Rcam'*inv(K)*q-tcam;
end
elseif quadric==2, %Parabolic Mirror
for i=1:n,    
    q = Q(:,i);
    v = Rcam'*inv(K)*q;
    Xhmir([1:2],i) =  [ v(1) ; v(2) ];
    Xhmir(3,i)=(Xhmir(1,i)^2+Xhmir(2,i)^2)/(2*a)-a/2;
    Xhmir(4,i)=1;%Homogeneous notation
end;    
end;    
