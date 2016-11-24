%%
% QuadricNLS minimizes Kd,R,t based on the cost function:
% ||C* - PQ*P||
%
% Input - Kd - Depth camera calibration matrix
%         Kr - RGB camera calibration matrix
%         Conic_RGB - Conics calculated from the RGB images of the spheres
%         R_hat - rotation matrix from R_R_D
%         t_hat - tranlstion matrix from R_t_D
%         U_depth - structure containing the points on the sphere
%         Switch - 0 for estimation of Kd,R,t
%                  1 for estimation of Kd
%
% Ouput - Kd_NLS - minimized Kd
%         R_NLS - minimized R
%         t_NLS - minimized t
%         RESIDUAL - residual from the non-linear least squares method
%
%%

function [Kd_NLS,R_NLS,t_NLS,RESIDUAL] = f_QuadricNLS(Kd, Kr, Conic_RGB, R_hat, t_hat, U_depth, Switch)


global U_depth_NLS Conic_RGB_NLS Kr_NLS Kd_NLS R_NLS t_NLS

U_depth_NLS = U_depth;
Conic_RGB_NLS = Conic_RGB;
Kr_NLS = Kr;
Kd_NLS = Kd;

R_NLS = R_hat;
t_NLS = t_hat;

fu_d = Kd(1,1);
fv_d = Kd(2,2);
u0_d = Kd(1,3);
v0_d = Kd(2,3);
skew_d = Kd(1,2);


tx = t_hat(1);
ty = t_hat(2);
tz = t_hat(3);

% Conversion of R into Euler angles
% Put E.angles in X0
[roll_est,pitch_est,yaw_est] = f_rotationMat2EulerAngles(R_hat);

% Minimization Kd

if Switch == 0,
    
    X = [roll_est;pitch_est;yaw_est;tx;ty;tz;fu_d;fv_d;u0_d;v0_d;skew_d];
    if ~isempty(strfind(version,'2011'))
    options = optimset('Algorithm','levenberg-marquardt','TolX',1e-6,'TolFun',1e-6,  ...
        'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');
    else 
options = optimset('LevenbergMarquardt','on','TolX',1e-6,'TolFun',1e-6,  ...
        'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');
    end
    
    lb = [];
    ub = [];
    
    [NLS,AA,RESIDUAL,BB,CC,DD,EE] = lsqnonlin(@f_minQuadric,X,lb,ub,options);
    
    Kd_NLS = [NLS(7) NLS(11)  NLS(9);
               0     NLS(8)   NLS(10);
               0      0        1];
%     Kd_NLS = Kd;
    R_NLS = rotoz(NLS(1))*rotoy(NLS(2))*rotox(NLS(3));
    
    t_NLS = [NLS(4);NLS(5);NLS(6)];
    
else
    
    X = [fu_d;fv_d;u0_d;v0_d;skew_d];
    
    if ~isempty(strfind(version,'2011'))
    options = optimset('Algorithm','levenberg-marquardt','TolX',1e-6,'TolFun',1e-6,  ...
        'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');
    else
    options = optimset('LevenbergMarquardt','on','TolX',1e-6,'TolFun',1e-6,  ...
        'Jacobian','off','MaxIter',400,'MaxFunEvals',700,'Display','iter');
    end
    
    lb = [];
    ub = [];
    
    [NLS,AA,RESIDUAL,BB,CC,DD,EE] = lsqnonlin(@f_minQuadric_KD,X,lb,ub,options);
    
    Kd_NLS = [NLS(1) NLS(5)  NLS(3);
        0     NLS(2)  NLS(4);
        0      0       1];
    
end

end



