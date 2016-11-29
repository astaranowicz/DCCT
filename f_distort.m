%%
%Calculates the Projection of the Sphere in the image from the estimated
%center of the ellipse
%
%Input -  XYZ      - normalized real-world coordinates
%         D_camera - the distortion vector for the camera
%                       form: [ k1 k2 k3 k4 k5];
%
%Output - XYZ      - distorted real-world coordinates
%%

function XYZ = f_distort(XYZ, D_camera)

assert(size(XYZ, 1) == 3, 'f_distort:DimensionMismatch', ...
    ['A parameter matrix of 3 rows is expected (', num2str(size(XYZ, 1)), ' received) for distort function.']);

XYZ(1:2,:) = XYZ(1:2,:) ./ XYZ(3,:); %Normalize

r2 = sum(XYZ(1:2,:).^2);
dr = 1.0 + D_camera(1) * r2 + D_camera(2) * r2.^2 + D_camera(5) * r2.^3;
dt = [2 * D_camera(3) * XYZ(1,:) .* XYZ(2,:) + D_camera(4) * (r2 + 2 * XYZ(1,:).^2);
      2 * D_camera(4) * XYZ(1,:) .* XYZ(2,:) + D_camera(3) * (r2 + 2 * XYZ(2,:).^2)];
XYZ(1:2,:) = (dr .* XYZ(1:2,:) + dt) .* XYZ(3,:);

end
