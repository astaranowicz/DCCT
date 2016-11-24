%%
% minQuadric is the cost function: ||C*-PQ*P||  to minimize Kd
%
%
% Input - Local - X - set of parameters to be minimized contains: Kd parameters
%         Global - U_depth_NLS - Depth camera, set of points on a sphere
%                  Conic_RGB_NLS - RGB camera, ellipse parametric parameters
%                  Kr_NLS - RGB camera calibration matrix
%                  R_NLS - rotation from R_R_D
%                  t_NLS - translation from R_t_D
%
% Output - F - residuals from the cost function
%
%%

function F = f_minQuadric_KD(X)

global U_depth_NLS Conic_RGB_NLS Kr_NLS R_NLS t_NLS

% Depth camera calibration matrix
Kd = [X(1) X(5)  X(3);
    0   X(2)  X(4);
    0    0     1];
%R_R_D
R = R_NLS;
%R_t_D
t = t_NLS;
% P- projection matrix from paper excluding the Kd
R_H_D = [R t];
%% Sphere fit to points to find the 3D center of the sphere
for i = 1:length(U_depth_NLS)
    %Converts the pixel center of the sphere to the 3D point
    X_depth(i).points = f_depth2XYZ(Kd,U_depth_NLS(i).points);
    % Linear Least Squares to find the model of the sphere
    M = f_sphereLinLS(X_depth(i).points(1:3,:));
    centerSphere_hat(i).center = M(1:3);
    radius_hat(i) = M(4);
    
    % Covariance factor which weight more the items closer than further
    % from the camera
    N(i)= norm(R*centerSphere_hat(i).center+t)^2;
end
%% Weighting for the cost function
N_max = max(N)+0.5;
N_min = min(N)-0.1;
% N_d = (N_max - N)./(N_min - N);
% std_d = std(N);
% W = 1- exp((-N_d.^2)/(std_d)^2);

N = (N_max - N)./(N_min - N);
W = 1- exp(-N.^2);

%keyboard
toleranceForImConic = 1;
%% Cost function using ||C*-PQ*P||
for i = 1:length(U_depth_NLS)
    % Depth camera - Conic
    Q = [eye(3) -centerSphere_hat(i).center;
        -centerSphere_hat(i).center'  centerSphere_hat(i).center'*centerSphere_hat(i).center-radius_hat(i)^2];
    
    temp_Conic = inv(Kr_NLS * R_H_D * inv(Q) * R_H_D' * Kr_NLS');
    
    conic_RGB = f_param2Conic_Ellipse(Conic_RGB_NLS(i).t(1),Conic_RGB_NLS(i).t(2),Conic_RGB_NLS(i).a,Conic_RGB_NLS(i).b,Conic_RGB_NLS(i).alpha);
    
    temp_Conic = temp_Conic/temp_Conic(3,3);
    conic_RGB = conic_RGB/conic_RGB(3,3);
    
    temp_f(i) = norm( conic_RGB- temp_Conic,'fro')*W(i);
    %Convert Conic to Parametric form
%     [x_t,y_t,a_t,b_t,alpha_t] = f_conic2Param(temp_Conic,toleranceForImConic);
%     %alpha_t is in radians
%     theta_t = alpha_t;%*180/pi;
%     
%     % RGB camera - Ellipse parametric parameters
%     a_rgb = Conic_RGB_NLS(i).a;
%     b_rgb = Conic_RGB_NLS(i).b;
%     theta_rgb = Conic_RGB_NLS(i).alpha;%*180/pi;
%     t_rgb = Conic_RGB_NLS(i).t;
    %Weighted
%     temp_f(i) = ( (abs(a_t-a_rgb)+abs(b_t-b_rgb))  +  (abs(x_t-t_rgb(1))+abs(y_t-t_rgb(2)))  +  abs(theta_t - theta_rgb)  )*W(i);
%     if i == 2
%         display(['Sphere',int2str(i)]);
%         figure(i)
%         f_plot_conicwparams(a_rgb,b_rgb, t_rgb(1),t_rgb(2), theta_rgb, 'g')
%         hold on
%         %         axis ij
%         %         axis equal
%         f_plot_conicwparams(a_t,b_t, x_t,y_t, theta_t, 'k')
%     end
end

F = temp_f;

end

