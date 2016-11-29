%%
% minQuadric is the cost function: ||C*-PQ*P||  to minimize Kd, Dd, R, t
%
%
% Input - Local - X - set of parameters to be minimized contains: Kd,Dd,R,t
%         Global - U_depth_NLS - Depth camera, set of points on a sphere
%                  Conic_RGB_NLS - RGB camera, ellipse parametric parameters
%                  Kr_NLS - RGB camera calibration matrix
%                  Dr_NLS - RGB camera distortion vector
%
% Output - F - residuals from the cost function
%
%%
function F = f_minQuadric(X)


global U_depth_NLS Conic_RGB_NLS Kr_NLS Dr_NLS Kd_NLS Dd_NLS Dd_NLS

% Depth camera calibration matrix
Kd = [X(7) X(11)  X(9);
    0   X(8)  X(10);
    0    0     1];
% Kd = Kd_NLS;
Dd = [X(12) X(13) X(14) X(15) X(16)];
% Dd = Dd_NLS;
%R_R_D
R = rotoz(X(1))*rotoy(X(2))*rotox(X(3));
%R_t_D
t = [X(4);X(5);X(6)];
% P- projection matrix from paper excluding the Kd and Dd
R_H_D = [R t];
%% Sphere fit to a set of points to find the 3D center of the sphere
for i = 1:length(U_depth_NLS)
    %Converts the pixel center of the sphere to the 3D point of the sphere
    X_depth(i).points = f_depth2XYZ(Kd,Dd,U_depth_NLS(i).points);
    % Linear Least Square to find the model of the sphere
    M = f_sphereLinLS(X_depth(i).points(1:3,:));
    centerSphere_hat(i).center = M(1:3);
    radius_hat(i) = M(4);
    
    % Covariance factor
    N(i)= norm(R*centerSphere_hat(i).center+t)^2;
end
%% Weighting for the cost function
N_max = max(N);
N_min = min(N);
N = (N_max - N)./(N_min - N);
W = 1- exp(-N.^2);

toleranceForImConic = 1e-6;
%% Cost function using ||C*-PQ*P||

try

for i = 1:length(U_depth_NLS)
    %Depth Camera Conic
    Q = [eye(3) -centerSphere_hat(i).center;
        -centerSphere_hat(i).center'  centerSphere_hat(i).center'*centerSphere_hat(i).center-radius_hat(i)^2];
    
    invQPt = inv(Q) * R_H_D' * Kr_NLS';
    invQPt(1:3,:) = f_undistort(invQPt(1:3,:) ./ invQPt(4,:), Dr_NLS) .* invQPt(4,:);

    temp_Conic = inv(Kr_NLS * R_H_D * invQPt);
    [x_t,y_t,a_t,b_t,alpha_t] = f_conic2Param(temp_Conic,toleranceForImConic);
    %alpha_t is in radians
    theta_t = alpha_t;%*180/pi;
    
    % RGB camera - Ellipse parametric parameters
    a_rgb = Conic_RGB_NLS(i).a;
    b_rgb = Conic_RGB_NLS(i).b;
    theta_rgb = Conic_RGB_NLS(i).alpha;%*180/pi;
    t_rgb = Conic_RGB_NLS(i).t;
    %Weighted
    temp_f(i) = (abs(a_t-a_rgb)+abs(b_t-b_rgb)+abs(x_t-t_rgb(1))+abs(y_t-t_rgb(2))+abs(theta_t - theta_rgb)) * W(i);
    
%     if i == 1
%         display(['Sphere',int2str(i)]);
%         figure(i)
%         f_plot_conicwparams(a_rgb,b_rgb, t_rgb(1),t_rgb(2), theta_rgb, 'g')
%         hold on
%         %         axis ij
%         %         axis equal
%         f_plot_conicwparams(a_t,b_t, x_t,y_t, theta_t, 'k')
% %     end
    
end

catch
	%TODO: Can you fix that?
	fprintf(2,'Optmization failed due to complex numbers in the conic parameters.\n')
	fprintf(2,'This is caused by the f_undistort function in this file\n')
	fprintf(2,'Terminating ...\n\n')
	return;
end

F = temp_f;

end


