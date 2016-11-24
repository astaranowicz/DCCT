function zx = imconicisect(p, w, t, R)
% Intersections of conic with window bounds.
%
% Input arguments:
% p:
%    the parameter vector p = [a b c d f g]
% w:
%    the window w = [x_min x_max y_min y_max]
% t:
%    the translation applied to the points
% R:
%    the rotation applied to the points
%
% zx:
%    points of intersection with window bounds in the standard coordinate
%    system (no rotation or translation applied)
%
% See also: imconic

% Copyright 2010 Levente Hunyadi

validateattributes(p, {'numeric'}, {'nonempty','real','vector'});
assert(length(p) == 6, 'imconicisect:DimensionMismatch', ...
    'A parameter vector of length 6 expected for a to g in a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0.');
if nargin < 2 || isempty(w)
    w = [];
else
    validateattributes(w, {'numeric'}, {'nonempty','real','vector'});
    assert(numel(w) == 4, 'imconicisect:DimensionMismatch', ...
        'A parameter vector of length 4 expected for x_min, x_max, y_min and y_max.');
end
if nargin < 4
    R = [];
    if nargin < 3
        t = [];
    end
end

if ~isempty(w)
    x_min = w(1); x_max = w(2);
    y_min = w(3); y_max = w(4);

    % calculate intersections with window bounds
    zx = [ ...
        imconicinterx(p, x_min, y_min, y_max) ...
        imconicinterx(p, x_max, y_min, y_max) ...
        imconicintery(p, y_min, x_min, x_max) ...
        imconicintery(p, y_max, x_min, x_max) ];
    %plot(zx(1,:),zx(2,:), 'ko');
    zx = invtransrot(zx,t,R);
else
    zx = zeros(2,0);
end

function z = imconicinterx(p, x, y_min, y_max)
% Intersection of conic section and straight line parallel to x-axis.

a = p(1); b = p(2); c = p(3); d = p(4); f = p(5); g = p(6);
eq = [c, b*x+f, a*x^2+d*x+g];  % equation for free variable y
y = quadeqrange(eq, y_min, y_max);
z = [ repmat(x, 1, numel(y)) ; y ];

function z = imconicintery(p, y, x_min, x_max)
% Intersection of conic section and straight line parallel to y-axis.

a = p(1); b = p(2); c = p(3); d = p(4); f = p(5); g = p(6);
eq = [a, b*y+d, c*y^2+f*y+g];  % equation for free variable x
x = quadeqrange(eq, x_min, x_max);
z = [ x ; repmat(y, 1, numel(x)) ];

function y = quadeqrange(eq, y_min, y_max)
% Compute solutions of quadratic equation in given range.
%
% Input arguments:
% eq:
%    the coefficients quadratic equations in descending order of powers
% y_min, y_max:
%    the range lower and upper bound

validateattributes(eq, {'numeric'}, {'nonempty','real','vector'});
validateattributes(y_min, {'numeric'}, {'real','scalar'});
validateattributes(y_max, {'numeric'}, {'real','scalar'});

a = eq(1); b = eq(2); c = eq(3);
if abs(a) < 1e-12  % a linear equation
    y = -c/b;
else  % a truely quadratic equation
    D = b^2-4*a*c;  % discriminant
    if D > 0
        sD = realsqrt(D);
        y = [(-b-sD) (-b+sD)]/(2*a);
    elseif D == 0
        y = (-b)/(2*a);
    else
        y = zeros(1,0);
    end
end
y = y(y >= y_min & y <= y_max);  % filter roots out of bounds

% y = roots(eq);
% y = reshape(y, 1, numel(y));  % row vector