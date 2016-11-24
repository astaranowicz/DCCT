function zrt = transrot(z, t, R)
% Apply translation and rotation to a set of points.
%
% See also: invtransrot

% Copyright 2010 Levente Hunyadi

if nargin > 2 && ~isempty(R)
    zr = R * z;
else
    zr = z;
end
if nargin > 1 && ~isempty(t)
    zrt = bsxfun(@plus, zr, t);
else
    zrt = zr;
end
