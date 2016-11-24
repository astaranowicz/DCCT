function D = imconicdiscr(p)
% Discriminant for conic section.
%
% See also: imconic

% Copyright 2010 Levente Hunyadi

a = p(1);
b = 0.5*p(2);
c = p(3);
D = det([ a b ; b c ]);
