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
%
% Output arguments:
% R:
%    rotation matrix
%
% See also: imconicrotate, imconictranslation

% Copyright 2010 Levente Hunyadi

a = p(1); b = p(2); c = p(3);
%assert(abs(p(4)) < 1e-12 && abs(p(5)) < 1e-12, 'imconicrotation:InvalidArgumentValue', 'Conic section should be centered in the origin.');

D = imconicdiscr(p);
if abs(D) > 1e-12
    if D >= 0  % ellipse
        ra = a; rc = c;
    else  % hyperbola
        ra = c; rc = a;
    end

    % compute angle enclosed with positive x-axis (angle of rotation)
    if b == 0  % conic is aligned with axes
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
        if a ~= c
            phi = 0.5*acot((a-c)/b);
            if ra > rc
                phi = 0.5*pi + phi;
            end
        else
            phi = -0.25*pi*sign(b);
        end
        sin_phi = sin(phi);
        cos_phi = cos(phi);
    end
else  % parabola
    if b == 0
        if p(1) ~= 0  % y = x^2 or y = -x^2
            phi = 0.5*pi;
            sin_phi = 1;
            cos_phi = 0;
        else  % if p(3) ~= 0  % x = y^2 or x = -y^2
            phi = 0;
            sin_phi = 0;
            cos_phi = 1;
        end
    else  % if b ~= 0, conic is rotated
        phi = 0.5*acot((a-c)/b);
        %phi = 0.5*atan(b/(a-c));
        sin_phi = sin(phi);
        cos_phi = cos(phi);
    end
end
R = [ cos_phi -sin_phi ; sin_phi cos_phi ];  % rotation matrix
