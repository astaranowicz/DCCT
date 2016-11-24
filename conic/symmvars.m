function [variable,type,sz] = symmvars(deftype, variable, varargin)
% Parse parameters for symbolic matrix creation.
% This function is used internally.
%
% See also: symm, symmd, symms

% Copyright 2010 Levente Hunyadi

validateattributes(deftype, {'char'}, {'nonempty','row'});
validatestring(deftype, {'positive','real','clear'});
error(nargchk(3, 256, nargin, 'struct'));
validateattributes(variable, {'char'}, {'nonempty','row'});
type = varargin{numel(varargin)};
if ischar(type)
    type = validatestring(type, {'real','positive','clear'});
    varargin(numel(varargin)) = [];  % remove last argument
else
    type = deftype;
end
switch numel(varargin)
    case 1  % symm([1 2 3]) or symm(4)
        dim = varargin{1};
        validateattributes(dim, {'numeric'}, {'positive','integer'});
        if isscalar(dim)  % symm(4)
            sz = [dim dim];
        else  % symm([1 2 3])
            validateattributes(dim, {'numeric'}, {'nonempty','vector'});
            sz = dim(:)';
        end
    otherwise  % symm(1,2,3,4)
        sz = zeros(1,numel(varargin));
        for k = 1 : numel(varargin)
            dim = varargin{k};
            validateattributes(dim, {'numeric'}, {'positive','integer','scalar'});
            sz(k) = dim;
        end
end

% only scalars, vector and matrices are supported but not 3d arrays
validateattributes(sz, {'numeric'}, {'size', [1 2]});
