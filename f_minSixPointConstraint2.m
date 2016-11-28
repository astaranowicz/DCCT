%%
%Non-Linear Least Squares for 6 point cost function to minimize only R, t
%
%Input - X - R - rotation from camera (R) to depth (D)
%            t - translation from camera (R) to depth (D)
%        Global -
%               U_depth_NLS - pixel points of the contour of the sphere in
%                       the depth
%               U_camera_spherecenter_NLS - pixel center of sphere in RGB
%                       camera
%               Switch_NLS - for weighted or non-weighted cost function
%               Kd_NLS - Depth camera calibraiton matrix
%               Kr_NLS - RGB camera calibration matrix
%
%Output - F -  vector of residual from both cost functions
%
%Cost Function for the 6point algorithm
%       0 = U_camera - Kr[R t]*X_depth
%
%  X = [fu_d; fv_d; u0_d; v0_d; tx; ty; tz; roll_est; pitch_est;yaw_est];
%%

function  F = f_minSixPointConstraint2(X)

global U_depth_NLS U_camera_spherecenter_NLS Switch_NLS Kd_NLS Kr_NLS Ellipse_C

% Depth camera Calibration matrix
Kd  = Kd_NLS;
% translation R_D
t = [X(1);X(2);X(3)];
% rotation R_D
R = rotoz(X(4))*rotoy(X(5))*rotox(X(6));
% RGB camera calibration matrix
% Kr = Kr_NLS;
Kr = [X(7) X(11) X(9);
       0   X(8) X(10);
       0    0   1];
Dr = [X(12) X(13) X(14) X(15) X(16)];
   
tolerance = 1e-8;
%% Sphere fit to set of points to find 3D center of sphere
for i = 1:length(U_depth_NLS)
    %Converts pixel center of sphere to 3D point
    X_depth(i).points= f_depth2XYZ(Kd,Dd,U_depth_NLS(i).points);
    %Sphere fitting for each sphere
    M = f_sphereLinLS(X_depth(i).points(1:3,:));
    centerSphere_hat(i).center = M(1:3);
    
    % Covariance factor which weight more the items closer than further
    % from the camera
    N(i)= norm(R*centerSphere_hat(i).center+t)^2;
    ProjectedCenter_camera(i).point = f_projectionSphere(Ellipse_C(i).t(1), Ellipse_C(i).t(2), Ellipse_C(i).a, Ellipse_C(i).b, Ellipse_C(i).alpha, Kr, Dr, tolerance);
end
%% Create Weights
N_max = max(N);
N_min = min(N);
N = (N_max - N)./(N_min - N);
W = 1- exp(-N.^2);

%% Cost Function for 6pnt Algorithm
%%Switch Between Weighted and Non-Weighted Cost Function
%%1 = non Weighted, 0 =  Weighted
if Switch_NLS
    for i = 1:length(U_depth_NLS)
        U_depth_temp = Kr*[R t] *[centerSphere_hat(i).center;1];
        U_depth_temp = U_depth_temp ./ U_depth_temp(3); 
        dist_6pnt(i) = norm([ProjectedCenter_camera(i).point;1] - U_depth_temp);
    end
else
    for i = 1:length(U_depth_NLS)
        U_depth_temp = Kr*[R t] *[centerSphere_hat(i).center;1];
        U_depth_temp = U_depth_temp ./ U_depth_temp(3);
        dist_6pnt(i) = norm([ProjectedCenter_camera(i).point;1] - U_depth_temp)*W(i);
        
    end
end

F = dist_6pnt;

end

