%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox v.1.3 (EGT)  %%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [lambda1,lambda2] = f_lambda(a,b,quadric,v)
%
%  f_lambda : Scale factors defining the intersection between a
%             mirror and a line going  passing through the origin.
%
%  Syntax:
%  ------
%  a,b = hyperbolic mirror parameters (scalar values).
%        HYPERBOLOID: (z+e)^2/a^2-(x^2+y^2)/b^2=1 where e=sqrt(a^2+b^2);
%        PARABOLOID:  z+a/2= (x^2+y^2)/(2a).
%  quadric = 1 for hyperbolic mirror;
%            2 for parabolic mirror.
%  v = column vector (3by1) with the coordinates, expressed in the mirror 
%      frame, of the point through which passed the line.
%
%  Descr: 
%  -----  This function returns the scale factors which define the
%         intersections between a two sheet hyperboloid of revolution 
%         with equation in the mirror frame
%              (y+e)^2/a^2-(x^2+z^2)/b^2=1 where e=sqrt(a^2+b^2) 
%         and a line passing from the origin of the mirror frame and 
%         from an other point defined as parameter (v).
%         The intersections are defined as lambda1*v and lambda2*v.
%
% Gian Luca Mariottini
%       November, 2005
%

function [lambda1,lambda2] = f_lambda(a,b,quadric,v);

if quadric==1,
    e = sqrt(a^2+b^2);
    if ((b^2*v(3)^2-a^2*v(1)^2-a^2*v(2)^2)==0),
        lambda1 = -(b^2)/(2*e*v(3));
        lambda2 = lambda1;
    else
        lambda1 = (b^2*(-e*v(3)+a*norm(v)))/(b^2*v(3)^2-a^2*v(1)^2-a^2*v(2)^2);
        lambda2 = (b^2*(-e*v(3)-a*norm(v)))/(b^2*v(3)^2-a^2*v(1)^2-a^2*v(2)^2);
    end;
elseif quadric==2,
    lambda1=a*(v(3) + norm(v))/(v(1)^2+v(2)^2);
    lambda2=a*(v(3) - norm(v))/(v(1)^2+v(2)^2);%Two solutions, only the first one for parab.mirrors.
end;    