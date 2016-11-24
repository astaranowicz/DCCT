function S = symm(varargin)
% Create a matrix of symbolic variables.
%
% See also: sym, symmd, symms
%
% Examples:
% A = symm('A',6,7)
%    creates a symbolic matrix
% M = symm('M',[6,7])
%    creates a symbolic matrix
% b = symm('b',7,1)
%    creates a symbolic vector
% A = symm('A',6,7,'real')
%    creates a symbolic matrix with real entries
% b = symm('b',7,1,'real')
%    creates a symbolic vector with real entries
% p = symm('b',7,1,'positive')
%    creates a symbolic vector with positive entries

% Copyright 2010 Levente Hunyadi

[variable,type,sz] = symmvars('clear',varargin{:});
S = sym;
S(sz(1),sz(2)) = 0;  % expand matrix to desired dimensions

% populate matrix with entries
if any(sz == 1)  % a vector, not a matrix
    if all(sz < 10)
        for k = 1 : max(sz)
            S(k) = sym([variable num2str(k, '%d')], type);
        end
    else
        for k = 1 : max(sz)
            S(k) = sym([variable '_' num2str(k, '%d')], type);
        end
    end
else
    if all(sz < 10)
        for j = 1 : sz(2)
            for i = 1 : sz(1)
                S(i,j) = sym([variable num2str(i, '%d') num2str(j, '%d')], type);
            end
        end
    else
        for j = 1 : sz(2)
            for i = 1 : sz(1)
                S(i,j) = sym([variable '_' num2str(i, '%d') '_' num2str(j, '%d')], type);
            end
        end
    end
end