%%
%Non-Linear Least Squares for 6 point cost function
%
%Input - Kd - calibration matrix for the depth camera
%        Dd - distortion vector for the depth camera
%        R - rotation from camera (R) to depth (D)
%        t - translation from camera (R) to depth (D)
%        Kr - calibration matrix for the RGB camera
%        Dr - distortion vector for the RGB camera
%        U_depth - pixel points of the contour of the sphere in the depth
%        U_camera - pixel points of the contour of the sphere in the camera
%        Switch - 0 for weighted cost function
%                 1 for non-weighted cost function
%        Switch2 - 1 for minimizing Kr,Dr,Kd,Dd,R,t
%                  0 for minimizing Kr,Dr,R,t
%
%
%Output - Kd_est  - estimated Depth camera calibration matrix
%         Dd_est - estimated Depth camera distortion vector
%         R_est - estimated R - rotation matrix
%         t_est - estimated t - translation
%         Kr_est - estimated RGB camera calibration matrix
%         Dr_est - estimated RGB camera distortion vector
%         RESIDUAL - the residual from the nonlinear least squares function
%%

function  [Kd_est,Dd_est,R_est,t_est,Kr_est,Dr_est,RESIDUAL] = f_sixPointConstraintNonLinLS(Kd, Dd, R, t, Kr, Dr, U_depth, U_camera_spherecenter,ellipse_camera, Switch, Switch2)

global U_depth_NLS U_camera_spherecenter_NLS Switch_NLS  Kd_NLS Dd_NLS Kr_NLS Dr_NLS Ellipse_C

Switch_NLS = Switch;
Ellipse_C = ellipse_camera;
U_depth_NLS = U_depth;
Kd_NLS = Kd;
Dd_NLS = Dd;
Kr_NLS = Kr;
Dr_NLS = Dr;
U_camera_spherecenter_NLS = U_camera_spherecenter;

fu_d = Kd(1,1);
fv_d = Kd(2,2);
u0_d = Kd(1,3);
v0_d = Kd(2,3);
skew_d = Kd(1,2);

k1_d = Dd(1);
k2_d = Dd(2);
k3_d = Dd(3);
k4_d = Dd(4);
k5_d = Dd(5);

fu_r = Kr(1,1);
fv_r = Kr(2,2);
u0_r = Kr(1,3);
v0_r = Kr(2,3);
skew_r = Kr(1,2);

k1_r = Dr(1);
k2_r = Dr(2);
k3_r = Dr(3);
k4_r = Dr(4);
k5_r = Dr(5);

tx = t(1);
ty = t(2);
tz = t(3);

% Conversion of R into Euler angles
% Put E.angles in X0
[roll_est,pitch_est,yaw_est] = f_rotationMat2EulerAngles(R);

%Minimization R,t,Kd
if Switch2 == 1
    X = [fu_d; fv_d; u0_d; v0_d; skew_d; k1_d;k2_d;k3_d;k4_d;k5_d; tx; ty; tz; roll_est; pitch_est;yaw_est; fu_r; fv_r; u0_r; v0_r;skew_r; k1_r;k2_r;k3_r;k4_r;k5_r];
    %     X = [fu_d; fv_d; u0_d; v0_d; skew_d; k1_d;k2_d;k3_d;k4_d;k5_d ;tx; ty; tz; roll_est; pitch_est;yaw_est];
	options = optimset('Algorithm','levenberg-marquardt','TolX',1e-6,'TolFun',1e-6,  ...
		'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');    
    
    lb = [];
    ub = [];
    
    [NLS,AA,RESIDUAL,BB,CC,DD,EE] = lsqnonlin(@f_minSixPointConstraint,X,lb,ub,options);
    
    Kd_est = [NLS(1)  NLS(5) NLS(3);
        0     NLS(2) NLS(4);
        0      0      1];
        
    Dd_est = [NLS(6) NLS(7) NLS(8) NLS(9) NLS(10)];
    
    t_est = [NLS(11);NLS(12);NLS(13)];
    R_est = rotoz(NLS(14))*rotoy(NLS(15))*rotox(NLS(16));
    
    Kr_est = [NLS(17)   NLS(21) NLS(19);
        0       NLS(18) NLS(20);
        0        0       1];
        
    Dr_est = [NLS(22) NLS(23) NLS(24) NLS(25) NLS(26)];
    
    %     Kr_est = Kr;
else
    
    X = [tx; ty; tz; roll_est; pitch_est;yaw_est;fu_r; fv_r; u0_r; v0_r;skew_r; k1_r;k2_r;k3_r;k4_r;k5_r];
	options = optimset('Algorithm','levenberg-marquardt','TolX',1e-6,'TolFun',1e-6,  ...
		'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');

    
    lb = [];
    ub = [];
    
    [NLS,AA,RESIDUAL,BB,CC,DD,EE] = lsqnonlin(@f_minSixPointConstraint2,X,lb,ub,options);
    
    t_est = [NLS(1);NLS(2);NLS(3)];
    R_est = rotoz(NLS(4))*rotoy(NLS(5))*rotox(NLS(6));
    %     Kr_est = Kr;
    Kr_est = [NLS(7)  NLS(11) NLS(9);
        0     NLS(8) NLS(10);
        0      0      1];
    %     Dr_est = Dr;
    Dr_est = [NLS(12) NLS(13) NLS(14) NLS(15) NLS(16)];
    Kd_est = Kd;    
    Dd_est = Dd;
    
end

end

