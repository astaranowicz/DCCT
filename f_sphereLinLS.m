%%
%Fits a sphere to a set of points using Linear Least Squares
%Input - World points on the sphere
%Output- Estimated Center and radius
%Algebraic equation of a Sphere:
% r^2 = (x-x0)^2+(y-y0)^2+(z-z0)^2
%
%Input -  Points - 3D points of the sphere
%
%Output - M - center(x,y) and radius
%
%%

function M = f_sphereLinLS(Points)

% Stacks the B matrix for SVD
B_hat = [(Points(1,:).^2 + Points(2,:).^2 + Points(3,:).^2)'  -2*Points(1,:)'  -2*Points(2,:)' -2*Points(3,:)' -ones(size(Points,2),1)];

% B_hat'*B_hat is less memory intensive than B_hat alone
[AA,BB,V] = svd(B_hat' * B_hat);

% x_hat = [1, x0, y0, z0, r^2-(x0^2+y0^2+z0^2)]
x_hat = V(:,5);
x_hat = x_hat/x_hat(1);

r = sqrt(x_hat(5) + (x_hat(2)^2 + x_hat(3)^2 + x_hat(4)^2));
center = [x_hat(2); x_hat(3); x_hat(4)];
radius_hat = r;

M = [center;radius_hat];

end

