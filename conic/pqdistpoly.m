function pqdistpoly
% Polynomial in parameter t for distance between point y and foot point x.
% Yields a polynomial in terms of the symbolic variable t, which is at most
% 4th degree for 2D and at most 6th degree for 3D. Solving the polynomial
% for t provides the minimum distance between the point y and its foot
% point x based on y - x = t * grad(Q(x)) = t * (2*A*x + b) and hence
% x = inv(I + 2*t*A) * (y - t*b).
% This function has been used to pregenerate the formulas in pqdist.
%
% See also: pqdist

% Copyright 2010 Levente Hunyadi

dim = 3;
alpha = symm('alpha',dim,1,'real');
beta = symm('beta',dim,1,'real');
c = sym('c','real');
t = sym('t','positive');
D = symmd('d',dim,dim,'positive');
d = diag(D);
M = symmd('m',dim,dim,'positive');  % M --> I + 2*t*D
m = diag(M);  % m(k) --> 1 + 2*t*d(k)

eq = (alpha-t*beta)'/M*D/M*(alpha-t*beta) + beta'/M*(alpha-t*beta) + c;  % = 0
eq = prod(m.^2) * eq;
eq = simplify(eq);  % needs m(k) to realize that expression can be simplified
for k = 1 : dim     % substitute m(k)
    eq = subs(eq, m(k), 1+2*t*d(k), 0);
end
[coeff,terms] = coeffs(eq, t);
for k = 1 : numel(terms)
    disp(terms(k));
    disp(sym2matlab(coeff(k)));
end
