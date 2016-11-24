function [R,phi] = imconicrotation(p)
% Rotation matrix of a conic section given as an implicit equation.
% The rotation matrix R rotates an axis-aligned conic section to the
% position as defined by p.
% The inverse of the rotation matrix R eliminates the xy term from the
% conic section equation.
%
% Input arguments:
% p:
%    parameters of the conic section implicit equation p = [a b c d f g]
% AS:
% tol: tolerance for the check if a != c
%
% Output arguments:
% R:
%    rotation matrix
%
% See also: imconicrotate, imconictranslation

% Copyright 2010 Levente Hunyadi
tol = 1e-8;
p=p/p(6);

a = p(1); b = p(2); c = p(3);
%assert(abs(p(4)) < 1e-12 && abs(p(5)) < 1e-12, 'imconicrotation:InvalidArgumentValue', 'Conic section should be centered in the origin.');

if imconicdiscr(p) > 0  % ellipse
    ra = a; rc = c;
else  % hyperbola
    ra = c; rc = a;
end
% compute angle enclosed with positive x-axis (angle of rotation)
if abs(b) < tol  % conic is aligned with axes
    if ra > rc
        phi = 0.5*pi;
        sin_phi = 1;
        cos_phi = 0;
    else  % ra < rc
       phi = 0;
       sin_phi = 0;
       cos_phi = 1;
    end
else  % if b ~= 0, conic is rotated
%          if a ~= c
    if abs(a-c) > tol  %AS:  the tolerance is a better soln since a&c are estimated
        phi = 0.5*acot((a-c)/b);
        if ra > rc
            phi = 0.5*pi + phi;
        end
    else
        if sign(b) == 1 ,  %AS: To check the orientation of the ellipse Note*.
            phi = 0.25*pi;
%             display('imconicrotation: pi/4');
        elseif sign(b) == -1,
            phi = -0.25*pi;
%             display('imconicrotation: -pi/4');
        end
    end
    sin_phi = sin(phi);
    cos_phi = cos(phi);
end
R = [ cos_phi -sin_phi ; sin_phi cos_phi ];  % rotation matrix



%AS: Note*: , the above check is due to special case that cos(pi/4) & sin(pi/4) ==
%cos(-pi/4) & sin(-pi/4)  since the conic parameters is determined by
%a, b, and c.
%when its pi/4, b = pos
%        -pi/4, b = neg

