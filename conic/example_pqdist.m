function example_pqdist
% Sample code for distance of quadratic curve and point.
%
% See also: pqdist

% Copyright 2010 Levente Hunyadi

% ellipse
p = imconicrotate([9 0 16 90 -64 145], pi/6);

P = transpose([0,1 ; -3,-6 ; -7,-3 ; -7,0]);  % each P = ( Px,Py )
[dst,F] = pqdist(P, p); %#ok<ASGLU>

figure('Name', 'Distance of points from ellipse');
hold('all');
imconic(p);
for k = 1 : size(F,2)
    plot([P(1,k) F(1,k)],[P(2,k) F(2,k)],'r');
end
axis equal
hold('off');

% hyperbola
p = imconicrotate([1 0 -1 0 0 -1], 7*pi/8);

P = transpose([0,1 ; -3,-6 ; -7,-3 ; -7,0 ; 10,5 ; 10,-10 ; 0,15]);  % each P = ( Px,Py )
[dst,F] = pqdist(P, p); %#ok<ASGLU>

figure('Name', 'Distance of points from hyperbola');
hold('all');
imconic(p);
for k = 1 : size(F,2)
    plot([P(1,k) F(1,k)],[P(2,k) F(2,k)],'r');
end
axis equal
hold('off');
