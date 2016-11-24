function [dst,X] = pqdist(Y, a, b, c)
% Distance of point(s) and quadratic curve or quadric surface.
% The quadratic curve (ellipse, parabola, hyperbola, etc.) or quadric
% surface (ellipsoid, elliptic paraboloid, hyperbolic paraboloid,
% hyperboloid, cone, elliptic cylinder, hyperbolic cylinder, parabolic
% cylinder, etc.) is given with the general quadratic equation
% Q(x) = x'*A*x + b'*x + c = 0
% where a pseudo-MatLab notation has been used, N = 2 for plane and N = 3
% for space. The parameter x is an N-by-1 column vector. Those points x
% that satisfy Q(x) = 0 comprise the quadratic curve or quadric surface.
%
% Input arguments:
% Y:
%    a column vector whose distance to the curve or surface to compute, or
%    a matrix of such column vectors
% -- -- --
% A:
%    a symmetric N-by-N matrix (N = 2 or N = 3 not necessarily invertible)
% b:
%    b is an N-by-1 column vector
% c:
%    a scalar
% -- OR --
% p:
%    a 6-element (for N = 2) or 10-element (for N = 3) parameter vector for
%    components x^2 xy y^2 x y 1 and x^2 y^2 z^2 x*y x*z y*z x y z 1,
%    respectively.
% -- -- --
%
% Output arguments:
% dst:
%    the distance between the point(s) Y and the curve or surface
% X:
%    a column vector of the foot point of Y on the curve or surface, or
%    a matrix of such column vectors
%
% References:
% David Eberly, Distance from Point to a General Quadratic Curve or a
%    General Quadric Surface,
%    http://www.geometrictools.com/Documentation/DistancePointToQuadratic.pdf

% Copyright 2010 Levente Hunyadi

validateattributes(Y, {'numeric'}, {'2d'});
if isvector(Y)  % a single point
    validateattributes(Y, {'numeric'}, {'column'}, mfilename, 'Y');
end

if nargin > 2
    validateattributes(a, {'numeric'}, {'2d'}, mfilename, 'A');
    A = a;
    dim = size(a,1);
    validateattributes(dim, {'numeric'}, {'>=', 2, '<=', 3});
    validateattributes(a, {'numeric'}, {'2d','size',[dim,dim]}, mfilename, 'A');
    validateattributes(b, {'numeric'}, {'column','size',[dim,1]}, mfilename, 'b');
    validateattributes(c, {'numeric'}, {'scalar'}, mfilename, 'c');
else
    validateattributes(a, {'numeric'}, {'nonempty','vector'}, mfilename, 'p');
    p = a;
    switch numel(p)
        case 6  % quadratic curve, x^2 xy y^2 x y 1
            dim = 2;
            A = [p(1) 0.5*p(2) ; 0.5*p(2) p(3)];
            b = [p(4) ; p(5)];
            c = p(6);
        case 10  % quadric surface, x^2 y^2 z^2 x*y x*z y*z x y z 1
            dim = 3;
            A = [ p(1), 0.5*p(4), 0.5*p(5) ...
                ; 0.5*p(4), p(2), 0.5*p(6) ...
                ; 0.5*p(5), 0.5*p(6), p(3) ...
                ];
            b = [p(7) ; p(8) ; p(9)];
            c = p(10);
        otherwise
            error('pqdist:DimensionMismatch', 'Unsupported parameter vector dimension numel(p) = %d.', numel(p));
    end
end

if rcond(A) < 1e-8  % degenerate quadratic curve/surface: line/plane
    %dst = b' * Y + c;  % row vector of signed distances
    %X = Y - b * dst;
    %return;
end

[R,D] = eig(A);
d = diag(D);
alpha = R'*Y;
beta = R'*b;

n = size(Y,2);
X = zeros(size(Y));
dst = zeros(1,n);

switch dim
    case 2  
        for k = 1 : n
            [dst(k),X(:,k)] = pqdistroots(Y(:,k), A, b, pqdist2(d, alpha(:,k), beta, c));
        end
    case 3  
        for k = 1 : n
            [dst(k),X(:,k)] = pqdistroots(Y(:,k), A, b, pqdist3(d, alpha(:,k), beta, c));
        end
end

function poly = pqdist2(d, alpha, beta, c)
% Distance polynomial of a single 2D point w.r.t. a quadratic curve.

poly = ...  % quadratic curve, polynomial coefficients precomputed with pqdistpoly
    [ - 4*beta(1)^2*d(1)*d(2)^2 - 4*beta(2)^2*d(1)^2*d(2) + 16*c*d(1)^2*d(2)^2 ... % t^4
    ; 4*beta(1)^2*d(1)*d(2) - beta(2)*(4*beta(2)*d(1)^2 + 8*beta(2)*d(2)*d(1)) - beta(1)*(4*beta(1)*d(2)^2 + 8*beta(1)*d(1)*d(2)) + 4*beta(2)^2*d(1)*d(2) + 16*c*d(1)*d(2)^2 + 16*c*d(1)^2*d(2) ... % t^3
    ; alpha(1)*(4*beta(1)*d(2)^2 + 8*beta(1)*d(1)*d(2)) + alpha(2)*(4*beta(2)*d(1)^2 + 8*beta(2)*d(2)*d(1)) - beta(1)*(2*beta(1)*d(1) + 4*beta(1)*d(2)) - beta(2)*(4*beta(2)*d(1) + 2*beta(2)*d(2)) + beta(1)^2*d(1) + beta(2)^2*d(2) + 4*c*d(1)^2 + 4*c*d(2)^2 + 4*alpha(1)^2*d(1)*d(2)^2 + 4*alpha(2)^2*d(1)^2*d(2) + 16*c*d(1)*d(2) - 8*alpha(1)*beta(1)*d(1)*d(2) - 8*alpha(2)*beta(2)*d(1)*d(2) ... % t^2
    ; 4*c*d(1) + 4*c*d(2) + alpha(1)*(2*beta(1)*d(1) + 4*beta(1)*d(2)) + alpha(2)*(4*beta(2)*d(1) + 2*beta(2)*d(2)) - beta(1)^2 - beta(2)^2 - 2*alpha(1)*beta(1)*d(1) - 2*alpha(2)*beta(2)*d(2) + 4*alpha(1)^2*d(1)*d(2) + 4*alpha(2)^2*d(1)*d(2) ... % t
    ; d(1)*alpha(1)^2 + beta(1)*alpha(1) + d(2)*alpha(2)^2 + beta(2)*alpha(2) + c ... % 1
    ];

function poly = pqdist3(d, alpha, beta, c)
% Distance polynomial of a single 3D point w.r.t. a quadric surface.

poly = ...  % quadric surface
    [ - 16*beta(1)^2*d(1)*d(2)^2*d(3)^2 - 16*beta(2)^2*d(1)^2*d(2)*d(3)^2 - 16*beta(3)^2*d(1)^2*d(2)^2*d(3) + 64*c*d(1)^2*d(2)^2*d(3)^2 ...
    ; 4*d(3)^2*(16*c*d(1)^2*d(2) + 16*c*d(1)*d(2)^2) - beta(1)*(4*d(3)^2*(4*beta(1)*d(2)^2 + 8*beta(1)*d(1)*d(2)) + 32*beta(1)*d(1)*d(2)^2*d(3)) - beta(2)*(4*d(3)^2*(4*beta(2)*d(1)^2 + 8*beta(2)*d(2)*d(1)) + 32*beta(2)*d(1)^2*d(2)*d(3)) + beta(1)^2*(16*d(1)*d(2)^2*d(3) + 16*d(1)*d(2)*d(3)^2) + beta(2)^2*(16*d(2)*d(1)^2*d(3) + 16*d(2)*d(1)*d(3)^2) + beta(3)^2*(16*d(3)*d(1)^2*d(2) + 16*d(3)*d(1)*d(2)^2) - beta(3)*(2*d(3)*(16*beta(3)*d(1)^2*d(2) + 16*beta(3)*d(1)*d(2)^2) + 16*beta(3)*d(1)^2*d(2)^2) + 64*c*d(1)^2*d(2)^2*d(3) ...
    ; alpha(2)*(4*d(3)^2*(4*beta(2)*d(1)^2 + 8*beta(2)*d(2)*d(1)) + 32*beta(2)*d(1)^2*d(2)*d(3)) + alpha(1)*(4*d(3)^2*(4*beta(1)*d(2)^2 + 8*beta(1)*d(1)*d(2)) + 32*beta(1)*d(1)*d(2)^2*d(3)) - beta(2)*(4*d(3)*(4*beta(2)*d(1)^2 + 8*beta(2)*d(2)*d(1)) + 4*d(3)^2*(4*beta(2)*d(1) + 2*beta(2)*d(2)) + 8*beta(2)*d(1)^2*d(2)) - beta(1)*(4*d(3)*(4*beta(1)*d(2)^2 + 8*beta(1)*d(1)*d(2)) + 4*d(3)^2*(2*beta(1)*d(1) + 4*beta(1)*d(2)) + 8*beta(1)*d(1)*d(2)^2) + 4*d(3)^2*(4*c*d(1)^2 + 16*c*d(1)*d(2) + 4*c*d(2)^2) + beta(1)^2*(4*d(1)*d(2)^2 + 16*d(1)*d(2)*d(3) + 4*d(1)*d(3)^2) + beta(2)^2*(4*d(2)*d(1)^2 + 16*d(2)*d(1)*d(3) + 4*d(2)*d(3)^2) + beta(3)^2*(4*d(3)*d(1)^2 + 16*d(3)*d(1)*d(2) + 4*d(3)*d(2)^2) + alpha(3)*(2*d(3)*(16*beta(3)*d(1)^2*d(2) + 16*beta(3)*d(1)*d(2)^2) + 16*beta(3)*d(1)^2*d(2)^2) - beta(3)*(2*d(3)*(4*beta(3)*d(1)^2 + 16*beta(3)*d(1)*d(2) + 4*beta(3)*d(2)^2) + 16*beta(3)*d(1)*d(2)^2 + 16*beta(3)*d(1)^2*d(2)) + 4*d(3)*(16*c*d(1)^2*d(2) + 16*c*d(1)*d(2)^2) + 16*c*d(1)^2*d(2)^2 - 2*alpha(1)*beta(1)*(16*d(1)*d(2)^2*d(3) + 16*d(1)*d(2)*d(3)^2) - 2*alpha(2)*beta(2)*(16*d(2)*d(1)^2*d(3) + 16*d(2)*d(1)*d(3)^2) - 2*alpha(3)*beta(3)*(16*d(3)*d(1)^2*d(2) + 16*d(3)*d(1)*d(2)^2) + 16*alpha(1)^2*d(1)*d(2)^2*d(3)^2 + 16*alpha(2)^2*d(1)^2*d(2)*d(3)^2 + 16*alpha(3)^2*d(1)^2*d(2)^2*d(3) ...
    ; alpha(2)*(4*d(3)*(4*beta(2)*d(1)^2 + 8*beta(2)*d(2)*d(1)) + 4*d(3)^2*(4*beta(2)*d(1) + 2*beta(2)*d(2)) + 8*beta(2)*d(1)^2*d(2)) + alpha(1)*(4*d(3)*(4*beta(1)*d(2)^2 + 8*beta(1)*d(1)*d(2)) + 4*d(3)^2*(2*beta(1)*d(1) + 4*beta(1)*d(2)) + 8*beta(1)*d(1)*d(2)^2) - beta(3)*(2*d(3)*(4*beta(3)*d(1) + 4*beta(3)*d(2)) + 4*beta(3)*d(1)^2 + 4*beta(3)*d(2)^2 + 16*beta(3)*d(1)*d(2)) - beta(2)*(4*d(3)*(4*beta(2)*d(1) + 2*beta(2)*d(2)) + 4*beta(2)*d(1)^2 + 4*beta(2)*d(3)^2 + 8*beta(2)*d(1)*d(2)) - beta(1)*(4*d(3)*(2*beta(1)*d(1) + 4*beta(1)*d(2)) + 4*beta(1)*d(2)^2 + 4*beta(1)*d(3)^2 + 8*beta(1)*d(1)*d(2)) + alpha(1)^2*(16*d(1)*d(2)^2*d(3) + 16*d(1)*d(2)*d(3)^2) + alpha(2)^2*(16*d(2)*d(1)^2*d(3) + 16*d(2)*d(1)*d(3)^2) + alpha(3)^2*(16*d(3)*d(1)^2*d(2) + 16*d(3)*d(1)*d(2)^2) + alpha(3)*(2*d(3)*(4*beta(3)*d(1)^2 + 16*beta(3)*d(1)*d(2) + 4*beta(3)*d(2)^2) + 16*beta(3)*d(1)*d(2)^2 + 16*beta(3)*d(1)^2*d(2)) + 4*d(3)*(4*c*d(1)^2 + 16*c*d(1)*d(2) + 4*c*d(2)^2) + 4*d(3)^2*(4*c*d(1) + 4*c*d(2)) + beta(1)^2*(4*d(1)*d(2) + 4*d(1)*d(3)) + beta(2)^2*(4*d(1)*d(2) + 4*d(2)*d(3)) + beta(3)^2*(4*d(1)*d(3) + 4*d(2)*d(3)) - 2*alpha(1)*beta(1)*(4*d(1)*d(2)^2 + 16*d(1)*d(2)*d(3) + 4*d(1)*d(3)^2) - 2*alpha(2)*beta(2)*(4*d(2)*d(1)^2 + 16*d(2)*d(1)*d(3) + 4*d(2)*d(3)^2) - 2*alpha(3)*beta(3)*(4*d(3)*d(1)^2 + 16*d(3)*d(1)*d(2) + 4*d(3)*d(2)^2) + 16*c*d(1)*d(2)^2 + 16*c*d(1)^2*d(2) ...
    ; alpha(1)^2*(4*d(1)*d(2)^2 + 16*d(1)*d(2)*d(3) + 4*d(1)*d(3)^2) + alpha(2)^2*(4*d(2)*d(1)^2 + 16*d(2)*d(1)*d(3) + 4*d(2)*d(3)^2) + alpha(3)^2*(4*d(3)*d(1)^2 + 16*d(3)*d(1)*d(2) + 4*d(3)*d(2)^2) + alpha(3)*(2*d(3)*(4*beta(3)*d(1) + 4*beta(3)*d(2)) + 4*beta(3)*d(1)^2 + 4*beta(3)*d(2)^2 + 16*beta(3)*d(1)*d(2)) + alpha(2)*(4*d(3)*(4*beta(2)*d(1) + 2*beta(2)*d(2)) + 4*beta(2)*d(1)^2 + 4*beta(2)*d(3)^2 + 8*beta(2)*d(1)*d(2)) + alpha(1)*(4*d(3)*(2*beta(1)*d(1) + 4*beta(1)*d(2)) + 4*beta(1)*d(2)^2 + 4*beta(1)*d(3)^2 + 8*beta(1)*d(1)*d(2)) + 4*d(3)*(4*c*d(1) + 4*c*d(2)) + beta(1)^2*d(1) + beta(2)^2*d(2) + beta(3)^2*d(3) + 4*c*d(1)^2 + 4*c*d(2)^2 + 4*c*d(3)^2 - beta(1)*(2*beta(1)*d(1) + 4*beta(1)*d(2) + 4*beta(1)*d(3)) - beta(2)*(4*beta(2)*d(1) + 2*beta(2)*d(2) + 4*beta(2)*d(3)) - beta(3)*(4*beta(3)*d(1) + 4*beta(3)*d(2) + 2*beta(3)*d(3)) + 16*c*d(1)*d(2) - 2*alpha(1)*beta(1)*(4*d(1)*d(2) + 4*d(1)*d(3)) - 2*alpha(2)*beta(2)*(4*d(1)*d(2) + 4*d(2)*d(3)) - 2*alpha(3)*beta(3)*(4*d(1)*d(3) + 4*d(2)*d(3)) ...
    ; 4*c*d(1) + 4*c*d(2) + 4*c*d(3) + alpha(1)*(2*beta(1)*d(1) + 4*beta(1)*d(2) + 4*beta(1)*d(3)) + alpha(2)*(4*beta(2)*d(1) + 2*beta(2)*d(2) + 4*beta(2)*d(3)) + alpha(3)*(4*beta(3)*d(1) + 4*beta(3)*d(2) + 2*beta(3)*d(3)) - beta(1)^2 - beta(2)^2 - beta(3)^2 + alpha(1)^2*(4*d(1)*d(2) + 4*d(1)*d(3)) + alpha(2)^2*(4*d(1)*d(2) + 4*d(2)*d(3)) + alpha(3)^2*(4*d(1)*d(3) + 4*d(2)*d(3)) - 2*alpha(1)*beta(1)*d(1) - 2*alpha(2)*beta(2)*d(2) - 2*alpha(3)*beta(3)*d(3) ...
    ; d(1)*alpha(1)^2 + beta(1)*alpha(1) + d(2)*alpha(2)^2 + beta(2)*alpha(2) + d(3)*alpha(3)^2 + beta(3)*alpha(3) + c ...
    ];

function [mindst,x] = pqdistroots(y, A, b, poly)
% Distance and foot point of a point w.r.t. a quadratic curve or quadric.

t = roots(poly);
t(abs(imag(t)) > 1e-16) = [];          % remove complex roots
dim = numel(y);
f = zeros(dim, numel(t));              % candidate foot points
dst = zeros(1, numel(t));
for k = 1 : numel(t)                   % compute foot point for each solution of t, f = inv(I + 2*t*A) * (y - t*b)
    f(:,k) = (eye(dim,dim) + 2*t(k)*A) \ (y - t(k)*b);
    dst(k) = sum((y - f(:,k)).^2);     % squared distance
end
[minsqdst,ix] = min(dst);
mindst = sqrt(minsqdst);
x = f(:,ix);                           % foot point with smallest distance to y
