%%
% Calculates the ratio of r^2/c^2 for the sphere center equation
% using the Parameteric to Conic conversion
%
% Input - x0,y0 - center of the ellipse
%         a,b - major, minor axis of the ellipse
%         alpha  - the tilt of the ellipse
%         K_camera - the calibration matrix of the camera
%         D_camera - the distortion vector of the camera
%
%
% Output - r_c - ratio of r^2/c^2
%
%
% from "Camera Calibration using Spheres" Agrawal '08
% Parameteric
%   c = cos(g)    s = sin(g)
%   x(t) = h + c[a * cost(t)] - s[b*sin(t)]
%   y(t) = k + c[a * cost(t)] - c[b*sin(t)]
% Conic
%   Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
%%

function  gamma_r = f_findRatioR_C(x0,y0,a,b,alpha, K_camera, D_camera, tolerance)

%Converts the Parametric form to the Conic
Conic = f_param2Conic_Ellipse(x0,y0,a,b,alpha);
% Conic = Conic/Conic(1,1);
%Need inverse to get the (3,3) which matches the r/c ratio
CS_star = inv(Conic);
C_star11 = CS_star(1,1);
C_star12 = CS_star(1,2);
C_star13 = CS_star(1,3);
C_star22 = CS_star(2,2);
C_star23 = CS_star(2,3);
C_star33 = CS_star(3,3);

u0 = K_camera(1,3);
v0 = K_camera(2,3);
fu = K_camera(1,1);
fv = K_camera(2,2);
%TODO: Utilize D_camera
% When C_star13 == u0 or when C_star23 == v0
if (abs((C_star13/C_star33)-u0) < tolerance) && (abs((C_star23/C_star33)-v0) > tolerance)    
    gamma_r(1) = 1;
    gamma_r(2) = inv(1- ((C_star33*fu^2)/(C_star11 - (C_star13*u0))));
elseif (abs((C_star13/C_star33)-u0) > tolerance) && (abs((C_star23/C_star33)-v0) < tolerance)
    gamma_r(1) = inv(1- ((C_star33*fv^2)/(C_star22 - (C_star23*v0))));
    gamma_r(2) = 1;
elseif (abs((C_star13/C_star33)-u0) < tolerance) && (abs((C_star23/C_star33)-v0) < tolerance)
    % When C_star13 == u0 && when C_star23 == v0
    gamma_r(1) = 1;
    gamma_r(2) = 1;
else
    %Default
    %New Version/simplified
    gamma_r(1) = 1 - ((((C_star12-(C_star23*u0)-(C_star13*v0))/ C_star33)+ (u0*v0))/(((C_star13/C_star33)-u0)*((C_star23/C_star33)-v0)));
    gamma_r(2) = 1 - ((((C_star12-(C_star23*u0)-(C_star13*v0))/ C_star33)+ (u0*v0))/(((C_star13/C_star33)-u0)*((C_star23/C_star33)-v0)));
end

end

%Old Versions
%r_c = 1 - ((C_star12/C_star33) - u0_camera*v0_camera+(-Cx + u0_camera)*v0_camera ...
%                 + (-Cy + v0_camera)*u0_camera)/((C_star13/C_star33 - u0_camera)*(C_star23/C_star33 - v0_camera));

%     %Simplified version of above
%     r_c = 1 - ((C_star12/C_star33) + u0_camera*v0_camera - Cx*v0_camera ...
%                     -Cy*u0_camera)/((C_star13/C_star33 - u0_camera)*(C_star23/C_star33 - v0_camera));

