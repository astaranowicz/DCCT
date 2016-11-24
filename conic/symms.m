function S = symms(varargin)
% Create a symmetric matrix of symbolic variables.
%
% See also: symm
%
% Examples:
% A = symms('A',3)
%    creates a symmetric symbolic 3-by-3 matrix
% A = symms('A',3,'positive')
%    creates a symmetric symbolic 3-by-3 matrix with positive entries

% Copyright 2010 Levente Hunyadi

[variable,type,sz] = symmvars('real', varargin{:});
type = validatestring(type, {'positive','real'});  % symmetric matrices cannot be complex
assert(numel(unique(sz)) == 1, 'symms:DimensionMismatch', 'A symmetric matrix must be square.');
sz = sz(1);

% populate matrix with entries
if sz == 1  % a scalar
    S = sym(variable);
else
    S = sym;
    S(sz,sz) = 0;  % expand matrix to desired dimensions
    if sz < 10
        for j = 1 : sz
            for i = j : sz
                S(i,j) = sym([variable num2str(i, '%d') num2str(j, '%d')], type);
                if i ~= j
                    S(j,i) = S(i,j);
                end
            end
        end
    else
        for j = 1 : sz
            for i = j : sz
                S(i,j) = sym([variable '_' num2str(i, '%d') '_' num2str(j, '%d')], type);
                if i ~= j
                    S(j,i) = S(i,j);
                end
            end
        end
    end
end