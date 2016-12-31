%%
%Calculates the Projection of the Sphere in the image from the estimated
%center of the ellipse
%
%Input -  XYZ      - distorted real-world coordinates
%         D_camera - the distortion vector for the camera
%                       form: [ k1 k2 k3 k4 k5];
%
%Output - XYZ      - normalized real-world coordinates
%%

function XYZ = f_normalize(XYZ, D_camera);

assert(size(XYZ, 1) == 3, 'f_normalize:DimensionMismatch', ...
    ['A parameter matrix of 3 rows is expected (', num2str(size(XYZ, 1)), ' received) for undistort function.']);

XYZ(1:2,:) = XYZ(1:2,:) ./ (ones(2, 1) * XYZ(3,:)); %Normalize

XY = XYZ(1:2,:);
for i = 1:20 % Newton's Iteration
	r2 = sum(XY.^2);
	dr =  1 + D_camera(1) * r2 + D_camera(2) * r2.^2 + D_camera(5) * r2.^3;
	dt = [2 * D_camera(3) * XY(1,:) .* XY(2,:) + D_camera(4) * (r2 + 2 * XY(1,:).^2);
	      2 * D_camera(4) * XY(1,:) .* XY(2,:) + D_camera(3) * (r2 + 2 * XY(2,:).^2)];
	XY = (XYZ(1:2,:) - dt) ./ (ones(2, 1) * dr);
end

XYZ(1:2,:) = XY .* (ones(2, 1) * XYZ(3,:));

end
