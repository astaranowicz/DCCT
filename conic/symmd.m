function S = symmd(varargin)
% Create a diagonal matrix of symbolic variables.
%
% See also: symm
%
% Examples:
% D = symmd('D',3)
%    creates a symbolic diagonal 3-by-3 matrix
% D = symmd('D',5,6)
%    creates a symbolic diagonal 5-by-6 matrix
% D = symmd('D',[7,12])
%    creates a symbolic diagonal 7-by-12 matrix

% Copyright 2010 Levente Hunyadi

[variable,type,sz] = symmvars('clear',varargin{:});
S = sym;
S(sz(1),sz(2)) = 0;  % expand matrix to desired dimensions

% populate matrix with entries
if min(sz) < 10
    for k = 1 : min(sz)
        S(k,k) = sym([variable num2str(k, '%d')], type);
    end
else
    for k = 1 : min(sz)
        S(k,k) = sym([variable '_' num2str(k, '%d')], type);
    end
end