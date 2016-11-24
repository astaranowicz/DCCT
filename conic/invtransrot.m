function ztr = invtransrot(z, t, R)
% Apply translation and inverse rotation to a set of points.
%
% See also: transrot

% Copyright 2010 Levente Hunyadi

if nargin > 1 && ~isempty(t)
    zt = bsxfun(@minus, z, t);
else
    zt = z;
end
if nargin > 2 && ~isempty(R)
    ztr = R \ zt;
else
    ztr = zt;
end
