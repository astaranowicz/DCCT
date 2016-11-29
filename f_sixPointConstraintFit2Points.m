%%
%Minimizing the Kd - calibration matrix of the depth camera
%               Dd - distortion vector for the depth camera
%               R - rotation from camera (R) to depth (D)
%               t - translation from camera (R) to depth (D)
%Using the cost function of the 6point algorithm for the R,t and Kd
%
%Input - Kd - calibration matrix for the depth camera
%        Dd - distortion vector for the depth camera
%        R - rotation from camera (R) to depth (D)
%        t - translation from camera (R) to depth (D)
%        Kr - calibration matrix for the RGB camera
%        Dr - distortion vector for the RGB camera
%        U_depth - pixel points of the contour of the sphere in the depth
%        U_camera - pixel points of the projected sphere center
%        Switch - 0 for weighted cost function
%                 1 for non-weighted cost function
%        Switch2 - 1 for minimizing Kr,Dr,Kd,Dt,R,t
%                  0 for minimizing Kr,Dr,R,t
%
%
%Output - bestKd - Depth camera calibration matrix
%         bestDd - Depth camera distortion vector
%         bestR - rotation from camera to depth
%         bestt - translation from camera to depth
%         bestKr - RGB camera calibraiton matrix
%         bestDr - RGB camera distortion vector
%         bestRes - residual from the consesus
%         bestinliers - set of inliers calculated from the consesus
%
%X = [fu_d; fv_d; u0_d; v0_d; tx; ty; tz; roll_est; pitch_est;yaw_est];
%%

function  [bestKd,bestDd,bestR,bestt,bestKr,bestDr,bestRes,bestinliers] = f_sixPointConstraintFit2Points(Kd, Dd, R, t, Kr, Dr, U_depth, U_camera_spherecenter,ellipse_camera, Switch, Switch2)

%% NLS 6Pnt using Kd on all inliers
[Kd_est,Dd_est,R_est,t_est,Kr_est,Dr_est,RESIDUAL] = ...
    f_sixPointConstraintNonLinLS(Kd, Dd, R, t, Kr, Dr, U_depth, U_camera_spherecenter,ellipse_camera, Switch, Switch2);

bestinliers = 1:length(U_depth);
bestKd = Kd_est;
bestDd = Dd_est;
bestR = R_est;
bestt = t_est;
bestKr = Kr_est;
bestDr = Dr_est;
bestRes = RESIDUAL;

end
