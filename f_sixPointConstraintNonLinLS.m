%%
%Non-Linear Least Squares for 6 point cost function
%
%Input - Kd - calibration matrix for the depth camera
%        R - rotation from camera (R) to depth (D)
%        t - translation from camera (R) to depth (D)
%        Kr - calibration matrix for the RGB camera
%        U_depth - pixel points of the contour of the sphere in the depth
%        U_camera - pixel points of the contour of the sphere in the camera
%        Switch - 0 for weighted cost function
%                 1 for non-weighted cost function
%        Switch2 - 1 for minimizing Kd,R,t
%                  0 for minimizing R,t
%
%
%Output - Kd_est  - estimated Depth camera calibration matrix
%         R_est - estimated R - rotation matrix
%         t_est - estimated t - translation
%         Kr_est - estimated RGB camera calibration matrix
%         RESIDUAL - the residual from the nonlinear least squares function
%%

function  [Kd_est,R_est,t_est,Kr_est,RESIDUAL] = f_sixPointConstraintNonLinLS(Kd, R, t, Kr, U_depth, U_camera_spherecenter,ellipse_camera, Switch, Switch2)

global U_depth_NLS U_camera_spherecenter_NLS Switch_NLS  Kd_NLS Kr_NLS Ellipse_C

Switch_NLS = Switch;
Ellipse_C = ellipse_camera;
U_depth_NLS = U_depth;
Kd_NLS = Kd;
Kr_NLS = Kr;
U_camera_spherecenter_NLS = U_camera_spherecenter;

fu_d = Kd(1,1);
fv_d = Kd(2,2);
u0_d = Kd(1,3);
v0_d = Kd(2,3);
skew_d = Kd(1,2);

fu_r = Kr(1,1);
fv_r = Kr(2,2);
u0_r = Kr(1,3);
v0_r = Kr(2,3);
skew_r = Kr(1,2);

tx = t(1);
ty = t(2);
tz = t(3);

% Conversion of R into Euler angles
% Put E.angles in X0
[roll_est,pitch_est,yaw_est] = f_rotationMat2EulerAngles(R);

%Minimization R,t,Kd
if Switch2 == 1
    X = [fu_d; fv_d; u0_d; v0_d; skew_d;tx; ty; tz; roll_est; pitch_est;yaw_est; fu_r; fv_r; u0_r; v0_r;skew_r];
    %     X = [fu_d; fv_d; u0_d; v0_d; skew_d;tx; ty; tz; roll_est; pitch_est;yaw_est];
    
    if ~isempty(strfind(version,'2011'))
        options = optimset('Algorithm','levenberg-marquardt','TolX',1e-6,'TolFun',1e-6,  ...
            'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');
    else
        options = optimset('LevenbergMarquardt','on','TolX',1e-6,'TolFun',1e-6,  ...
            'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');
    end
    
    lb = [];
    ub = [];
    
    [NLS,AA,RESIDUAL,BB,CC,DD,EE] = lsqnonlin(@f_minSixPointConstraint,X,lb,ub,options);
    
    Kd_est = [NLS(1)  NLS(5) NLS(3);
        0     NLS(2) NLS(4);
        0      0      1];
    
    t_est = [NLS(6);NLS(7);NLS(8)];
    R_est = rotoz(NLS(9))*rotoy(NLS(10))*rotox(NLS(11));
    
    Kr_est = [NLS(12)   NLS(16) NLS(14);
        0       NLS(13) NLS(15);
        0        0       1];
    
    %     Kr_est = Kr;
else
    
    X = [tx; ty; tz; roll_est; pitch_est;yaw_est;fu_r; fv_r; u0_r; v0_r;skew_r];
    
    if ~isempty(strfind(version,'2011'))
        options = optimset('Algorithm','levenberg-marquardt','TolX',1e-6,'TolFun',1e-6,  ...
            'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');
    else
        options = optimset('LevenbergMarquardt','on','TolX',1e-6,'TolFun',1e-6,  ...
            'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');
    end
    
    lb = [];
    ub = [];
    
    [NLS,AA,RESIDUAL,BB,CC,DD,EE] = lsqnonlin(@f_minSixPointConstraint2,X,lb,ub,options);
    
    t_est = [NLS(1);NLS(2);NLS(3)];
    R_est = rotoz(NLS(4))*rotoy(NLS(5))*rotox(NLS(6));
    %     Kr_est = Kr;
    Kr_est = [NLS(7)  NLS(11) NLS(9);
        0     NLS(8) NLS(10);
        0      0      1];
    Kd_est = Kd;
    
end

end

