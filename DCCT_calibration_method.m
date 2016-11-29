%%
%
%%

function DCCT_calibration_method(Switch,varargin)
global DCCT_variables

pause(0.05)
if Switch == 1
    load Data/DepthSpherePoints
    load Data/RGBPixelcenters
elseif Switch == 2
    load Depth_userInput
    load RGB_userInput
else
    avg_depthPixelPoints = DCCT_variables.avg_depthPixelPoints;
    projected_sphere_center_RGB = DCCT_variables.projected_sphere_center_RGB;
    U_depth_clicked = DCCT_variables.U_depth_clicked;
    Ellipse_RGB = DCCT_variables.Ellipse_RGB;
end

%% 1) 6pnt Algorithm using average centers spheres -> (Kd, R, t)
display('Least Squares Phase')
guessInliers = DCCT_variables.setOfSpheres;
[Kd_hat, R_hat, t_hat, res_6pnt_avgcenter, Inliers_6pnt] = ...
    f_DepthCal_points([projected_sphere_center_RGB(guessInliers).points], avg_depthPixelPoints(:,guessInliers), DCCT_variables.Kr, DCCT_variables.threshold_Res6pnt);


%% 2) SphereFit (w/ est Kd) --> inliers_pnts(j)
for j = DCCT_variables.setOfSpheres
    X_depth(j).points = f_depth2XYZ(Kd_hat, U_depth_clicked(guessInliers(j)).points);
    [centerSphere_hat(j).center, radiusSphere_hat(j), residualSphere_hat(j).res, inliers_est(j).points, outliers_est(j).points, ind_inlsphere(j).index] = ...
        f_sphereFit_points2Sphere(X_depth(j).points(1:3,:), DCCT_variables.threshold_3D_Depth);
    % Depth-map pnts without outliers
    U_depth_inl(guessInliers(j)).points =  U_depth_clicked(guessInliers(j)).points(:,ind_inlsphere(j).index);
end


%% NLS
display('Nonlinear Minimization Phase')
display('6pnt Nonlinear Least Squares')
Switch = 0; % 0 for weighted and 1 for non-weighted
Switch2 = 1; % 0 to minimize Kr,R,t and 1 to minimise Kd,R,t
[Kd_C,R_C,t_C,Kr_C,RESNLS,inliers_6pnt_NLS] = f_sixPointConstraintFit2Points(Kd_hat, R_hat, t_hat,DCCT_variables.Kr, U_depth_inl(guessInliers(Inliers_6pnt)),...
    projected_sphere_center_RGB(guessInliers(Inliers_6pnt)),Ellipse_RGB(guessInliers(Inliers_6pnt)),Switch,Switch2);

display('Quadric Nonlinear Least Squares')
Switch = 1; % 0 for Kd,R,t and 1 for Kd
[Kd_C2,R_C2,t_C2,RES_NLS] = f_QuadricNLS(Kd_C, Kr_C, Ellipse_RGB(guessInliers(Inliers_6pnt)),R_C,t_C,...
    U_depth_inl(guessInliers(Inliers_6pnt)),Switch);

inliers_6pnt_NLS2=inliers_6pnt_NLS;
Kd_D2=Kd_C2;
R_D2=R_C2;
t_D2=t_C2;
Kr_D2=Kr_C;

%% Residual 6pnt After NLS minimization
display('Plots')
figure
hold on
title('Residual')
for i = 1:length(guessInliers(Inliers_6pnt))
    %Calculates the projected sphere center
    j =(guessInliers(Inliers_6pnt(i)));
    projected_sphere_center_inl_RGB_new(i).points = f_projectionSphere( Ellipse_RGB(j).t(1),  Ellipse_RGB(j).t(2),  Ellipse_RGB(j).a,  Ellipse_RGB(j).b,  Ellipse_RGB(j).alpha, Kr_D2,2);
    
    pnts(:,i) = f_residual6pntCostFunction(Kr_D2,Kd_D2,R_D2,t_D2,projected_sphere_center_inl_RGB_new(i), U_depth_inl(guessInliers(Inliers_6pnt(i))));
    plot(pnts(1,i),pnts(2,i),'rx');
    err(i) = norm(pnts(:,i));
end

meanErr = mean(err);
stdErr = std(err);

%% Displays the Projected Conic
Projected_Conic = f_Quadric2Conic(U_depth_inl,Kd_D2,Kr_D2,R_D2,t_D2);

for i = DCCT_variables.setOfSpheres
    selected_Conic = f_param2Conic_Ellipse(Ellipse_RGB(i).t(1),Ellipse_RGB(i).t(2),Ellipse_RGB(i).a,Ellipse_RGB(i).b,Ellipse_RGB(i).alpha);
    
    imageName = [DCCT_variables.dirExt,DCCT_variables.RGBImageName,num2str(i),DCCT_variables.fileExt];
    imageRead = imread(imageName);
    figure
    imshow(imageRead);
    hold on
    check =  find(i == guessInliers(Inliers_6pnt));
    if isempty(check) == 1
        title(['Sphere ',num2str(i),'is outlier'])
    else
        title(['Sphere ',num2str(i),'is inlier'])
    end
    f_drawConic(Projected_Conic(i).conic,1,'b')
    f_drawConic(selected_Conic,1,'g');
end
display(' ')
display(' ')
display(' ')
display('Depth sensor calibration matrix')
display(num2str(Kd_C2))
display('RGB sensor calibration matrix')
display(num2str(Kr_C))
display('Extrinsic calibration parameters: R')
display(num2str(R_C2))
display('Extrinsic calibration parameters: t [meters]')
display(num2str(t_C2))
display(['Pixel Reprojection Error: ',num2str(meanErr),' +(-) ',num2str(stdErr), '(mean and standard deviation)' ])
display(' ')
display(' ')
userSaveRes = input('Do you want to save the calibration results: y,n [n]:  ','s');
if ~isempty(strfind(userSaveRes,'y'))
    display('Saving');
    save('CalibrationResults','Kd_C2','Kr_C','R_C2','t_C2','err');
else
    display('Not saving');
end



end




