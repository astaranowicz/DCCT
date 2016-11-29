%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox  v1.3 (EGT)  %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  f_projectmir(a,b,quadric,Xmir) 
%
%  Usage:
%  -----
%      Projection on the mirror surface of a point Xmir in the mirror frame.
%
%  Sintax:
%  ------
%  a,b = mirror parameters:
%        HYPERBOLOID: (z+e)^2/a^2-(x^2+y^2)/b^2=1 where e=sqrt(a^2+b^2);
%        PARABOLOID:  z+a/2= (x^2+y^2)/(2a).
%  quadric = 1 for hyperbolic mirror;
%            2 for parabolic mirror.
%  Xmir = coordinates of point expressed in the mirror frame
%
%  Description: 
%  -----------
%         This function projects a point, whose coordinates are expressed 
%         in the mirror coordinate system, on the mirror surface. 
%         Also the projection on the mirror surface is expressed in the
%         mirror coordinates system.
%
%
% Author:
%     Gian Luca Mariottini
% Last update:
%     Dec, 2005

function [pmir] = f_projectmir(a,b,quadric,Xmir);

[lambda1,lambda2] = f_lambda(a,b,quadric,Xmir);
if quadric==1,
    if lambda1==lambda2 ,
        pmir = lambda1*Xmir;
    else
    if Xmir(3)<0,
        r = abs(Xmir(3))*(b/a);
        if (abs(Xmir(1))<r)&(abs(Xmir(2))<r),
            lambda=min([lambda1,lambda2]);
            pmir=lambda*Xmir;
        else 
            lambda=max([lambda1,lambda2]);
            pmir=lambda*Xmir;
        end
    else 
        lambda=min([lambda1,lambda2]);
        pmir=lambda*Xmir;
    end;
  end;
elseif quadric==2,
    lambda=lambda1;
    pmir=lambda*Xmir;
end