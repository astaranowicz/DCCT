function example_symm
% Sample code for matrices in symbolic variables
%
% See also: symm, symmd, symms

% Copyright 2010 Levente Hunyadi

dim = 3;
A = symms('A',dim,dim,'real');
b = symm('b',dim,1,'real');
c = sym('c','real');
x = symm('x',dim,1,'real');

Q = expand(x'*A*x + b'*x + c);

% Expanded terms for dim = 2 are as follows:
% A(1,1)*x(1)^2
% A(2,2)*x(2)^2
% 2*A(2,1)*x(1)*x(2)
% b(1)*x(1)
% b(2)*x(2)
% c

% Expanded terms for dim = 3 are as follows:
% A(1,1)*x(1)^2
% A(2,2)*x(2)^2
% A(3,3)*x(3)^2
% 2*A(2,1)*x(1)*x(2)
% 2*A(3,1)*x(1)*x(3)
% 2*A(3,2)*x(2)*x(3)
% b(1)*x(1)
% b(2)*x(2)
% b(3)*x(3)
% c

disp(sym2matlab(Q));
