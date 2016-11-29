%%
%Function is used to determine if the points are within a certain threshold
%on the estimated model.
%Residual for a sphere follows this equation:
% d_i = ||Xi - X0|| - r
%
%Input - M - (x,y,z) center and (r) radius
%        x - points of the sphere
%        t - threshold to check if the points lie on the sphere or not
%
%Output - inliers - the points on the sphere with a certain threshold
%         M - (x,y,z) center and (r) radius
%         outliers - the points not on the sphere
%%

function [inliers, M, outliers,indices] = f_sphereResidual(M,x,t)

inliers = [];
outliers = [];

x0 = M(1); % x center
y0 = M(2); % y center
z0 = M(3); % z center
radius_i = M(4);  % radius of sphere

x2 = x- [x0; y0; z0] * ones(1, size(x,2));

dist = f_roundn(sqrt(sum(x2 .* x2)) - radius_i,-10);

indices = dist.^2 <= t^2;
inliers = x(:,indices);
outliers = x(:,logical(1-indices));

end




