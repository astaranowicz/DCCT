%%
%Calculates the Projection of the Sphere in the image from the estimated
%center of the ellipse
%
%Input -  XYZ - distorted real-world coordinates
%         D_camera - the distortion vector for the camera
%                       form: [ k1 k2 k3 k4 k5];
%
%Output - XYZ - undistorted real-world coordinates
%%

function XYZ = f_undistort(XYZ, D_camera)

assert(size(XYZ, 1) >= 2 && size(XYZ, 1) <= 3, 'f_undistort:DimensionMismatch', ...
    'A parameter matrix of 2 or 3 rows is expected for undistort function.');

r2 = XYZ(1,:).^2 + XYZ(2,:).^2;
dr = 1 + D_camera(1) * r2 + D_camera(2) * r2.^2 + D_camera(5) * r2.^3;
dtx = 2 * D_camera(3) * XYZ(1,:) .* XYZ(2,:) + D_camera(4) * (r2 + 2 * XYZ(1,:).^2);
dty = 2 * D_camera(4) * XYZ(1,:) .* XYZ(2,:) + D_camera(3) * (r2 + 2 * XYZ(2,:).^2);
XYZ(1,:) = dr .* XYZ(1,:) + dtx;
XYZ(2,:) = dr .* XYZ(2,:) + dty;

end
