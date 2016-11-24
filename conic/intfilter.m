function [z,ix] = intfilter(x, bnd)
% Filters those elements of x that are outside a set of intervals.
%
% Input arguments:
% x:
%    a vector whose elements to check if they are in the set of intervals
% bnd:
%    a matrix whose first row is the lower bound and second row is the
%    upper bound for the set of intervals
%
% Output arguments:
% ix:
%    in index of the interval in which the respective elements of x have
%    been found

% Copyright 2010 Levente Hunyadi

validateattributes(x, {'numeric'}, {'real','nonempty','vector'});
validateattributes(bnd, {'numeric'}, {'real','2d'});
assert(size(bnd,1) == 2, 'isininterval:DimensionMismatch', ...
    'Interval bounds should be supplied as a 2-by-n matrix.');

isrow = size(x,2) > size(x,1);
x = x(:);  % column vector, x_i with i = 1..numel(x)

I = bsxfun(@minus, x, bnd(1,:)) >= 0 & bsxfun(@minus, bnd(2,:), x) >= 0;  % indicator matrix, 1 if x_i is in interval j
if nargout > 1
    ix = zeros(numel(x),1);
    for k = 1 : numel(x)
        ix(k) = max([ 0 find(I(k,:), 1) ]);  % index of first interval x_i is in (or 0 if none)
    end
    z = x(ix > 0);
    ix = ix(ix > 0);
    if isrow
        ix = ix.';
    end
else
    z = x(any(I,2));  % elements of x that are in any of the intervals
end
    
if isrow
    z = z.';
end