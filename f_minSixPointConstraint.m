%%
%Non-Linear Least Squares for 6point cost function to minimize Kd, R, and t
%
%Input - X - Kd - calibration matrix for the depth camera
%            R - rotation from camera (R) to depth (D)
%            t - translation from camera (R) to depth (D)
%        Global -
%               U_depth_NLS - pixel points of the contour of the sphere in the depth
%               U_camera_spherecenter_NLS - pixel center of sphere in RGB
%                       camera
%               Switch_NLS - for weighted/non-weighted cost function
%               Kr_NLS - RGB camera calibration matrix
%
%Output - F -  vector of residual from both cost functions
%
%Cost Function for the 6point algorithm
%       0 = U_camera - Kr[R t]*X_depth
%
%  X = [fu_d; fv_d; u0_d; v0_d; tx; ty; tz; roll_est; pitch_est;yaw_est];
%%

function  F = f_minSixPointConstraint(X)

global U_depth_NLS U_camera_spherecenter_NLS Switch_NLS Kr_NLS Ellipse_C

% Depth camera calibration matrix
Kd  = [X(1) X(5) X(3);
    0   X(2) X(4);
    0    0   1];
% R_t_D
t = [X(6);X(7);X(8)];
% R_R_D
R = rotoz(X(9))*rotoy(X(10))*rotox(X(11));
% RGB camera calibration matrix
% Kr = Kr_NLS;
Kr = [X(12) X(16) X(14);
    0   X(13) X(15);
    0    0   1];
%% Sphere fit to a set of points to find the 3D center of sphere
tolerance = 1e-8;
for i = 1:length(U_depth_NLS)
    %Converts pixel center of sphere to 3D point
    X_depth(i).points= f_depth2XYZ(Kd,U_depth_NLS(i).points);
    
    %Sphere fitting for each sphere
    M = f_sphereLinLS(X_depth(i).points(1:3,:));
    centerSphere_hat(i).center = M(1:3);
    
    % Covariance factor which weight more the items closer than further
    % from the camera
    N(i)= norm(R*centerSphere_hat(i).center+t)^2;
    ProjectedCenter_camera(i).point = f_projectionSphere(Ellipse_C(i).t(1), Ellipse_C(i).t(2), Ellipse_C(i).a, Ellipse_C(i).b, Ellipse_C(i).alpha, Kr,tolerance);
end
%% Create Weights
N_max = max(N)+0.5;
N_min = min(N)-0.1;
% N_d = (N_max - N)./(N_min - N);
% std_d = std(N);
% W = 1- exp((-N_d.^2)/(std_d)^2);

N = (N_max - N)./(N_min - N);
W = 1- exp(-N.^2);
%% Cost Function for 6pnt Algorithm
%%Switch Between Weighted and Non-Weighted Cost Function
%%1 = non Weighted, 0 =  Weighted
if Switch_NLS
    for i = 1:length(U_depth_NLS)
        U_depth_temp = Kr*[R t] *[centerSphere_hat(i).center;1];
        U_depth_temp = U_depth_temp ./ U_depth_temp(3);
%         dist_6pnt(i) = norm([U_camera_spherecenter_NLS(i).points;1] - U_depth_temp);
        dist_6pnt(i) = norm([ProjectedCenter_camera(i).point;1] - U_depth_temp);
    end
else
    for i = 1:length(U_depth_NLS)
        U_depth_temp = Kr*[R t] *[centerSphere_hat(i).center;1];
        U_depth_temp = U_depth_temp ./ U_depth_temp(3);
%         dist_6pnt(i) = norm([U_camera_spherecenter_NLS(i).points;1] - U_depth_temp)*W(i);
        dist_6pnt(i) = norm([ProjectedCenter_camera(i).point;1] - U_depth_temp)*W(i);
    end
end

F = dist_6pnt;

end
