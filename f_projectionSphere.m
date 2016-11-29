%%
%Calculates the Projection of the Sphere in the image from the estimated
%center of the ellipse
%
%Input - x0,y0 - center of the ellipse
%        a,b  - major, minor axis of the ellipse
%        alpha - tilt of the ellipse
%        K_camera - the calibration matrix for the camera
%                       form: [ fku 0  u0;
%                                0 fkv v0;
%                                0  0   1];
%        D_camera - the distortion vector for the camera
%                       form: [ k1 k2 k3 k4 k5];
%
%Output - sphere_center - (u,v) of the projected sphere center from the
%                               ellipse center
%%

function sphere_center = f_projectionSphere(x0, y0, a, b, alpha, K_camera, D_camera, tolerance)

%To find the ratio of R/C to be used to find the sphere projected center
gamma_r = f_findRatioR_C(x0, y0, a, b, alpha, K_camera, D_camera, tolerance);

%The center of the projected sphere in the camera frame
%TODO: Utilize D_camera
u = x0 * (1 - gamma_r(1)) + K_camera(1,3) * gamma_r(1);
v = y0 * (1 - gamma_r(2)) + K_camera(2,3) * gamma_r(2);

sphere_center = [u;v];
%  sphere_center(1) =  (1 - r_c) * (x0 - K_camera(1,3));
%  sphere_center(2) =  (1 - r_c) * (y0 - K_camera(2,3));

end
