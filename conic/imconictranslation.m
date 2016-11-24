function t = imconictranslation(p)
% Translation vector of a conic section given as an implicit equation.
% The translation vector t maps a (possibly rotated) conic section in the
% origin to the location as defined by p.
%
% Input arguments:
% p:
%    parameters of the conic section implicit equation p = [a b c d f g]
%
% Output arguments:
% t:
%    translation (column) vector
%
% See also: imconictranslate, imconicrotation

% References:
% "Classification and Principal Axis transformation of quadratic curves",
%    http://www.math4all.in/public_html/linear%20algebra/chapter11.3.html

% Copyright 2010 Levente Hunyadi

A_xx = p(1);
A_xy = 0.5*p(2);
A_yy = p(3);
B_x  = 0.5*p(4);
B_y  = 0.5*p(5);

D = imconicdiscr(p);
if abs(D) > 1e-12  % ellipse or hyperbola
    x0 = -det([ B_x A_xy ; B_y A_yy ]) / D;
    y0 = -det([ A_xx B_x ; A_xy B_y ]) / D;
else  % parabola
    % equation: a*x^2 + b*x*y + c*y^2 + p*x + q*y + d = (alpha*x + beta*y)^2 + p*x + q*y + d = 0
    [~,theta] = imconicrotation(p);
    p = imconicrotate(p, -theta);

    if abs(p(3)) > abs(p(1))
        % x = a*y^2 + b*y + c, p = [0 0 a 1 b c]
        p = p ./ p(4);
        a = p(3); b = p(5); c = p(6);
        x0 = (4*a*c - b^2) / (4*a);
        y0 = -b / (2*a);
    else
        % y = a*x^2 + b*x + c, p = [a 0 0 b 1 c]
        p = p ./ p(5);
        a = p(1); b = p(4); c = p(6);
        x0 = -b / (2*a);
        y0 = (4*a*c - b^2) / (4*a);
    end
end
t = [ x0 ; y0 ];  % translation vector

% use implicit equation to compute circle center
% x0 = -d/(2*a);
% y0 = -e/(2*a);
% t = [ x0 ; y0 ];
% use implicit equation to compute ellipse center
% x0 = (c*d-b*f)/(b^2-a*c);
% y0 = (a*f-b*d)/(b^2-a*c);
% t = [ x0 ; y0 ];  % translation vector
