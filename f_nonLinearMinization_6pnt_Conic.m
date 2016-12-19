%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C) 2016 Aaron Staranowicz and Gian Luca Mariottini
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%
%%

function Sol_NLS = f_nonLinearMinization_6pnt_Conic(Sol_LS, f_Opt)
% Extract parameters to optimize
[r0,p0,y0] = f_R2rpy(Sol_LS.R_R_D);
tx = Sol_LS.R_t_D(1); ty = Sol_LS.R_t_D(2); tz = Sol_LS.R_t_D(3);
D_fu=Sol_LS.K_D(1,1); D_fv=Sol_LS.K_D(2,2); D_s=0*Sol_LS.K_D(1,2); D_u0=Sol_LS.K_D(1,3); D_v0=Sol_LS.K_D(2,3);
R_fu=Sol_LS.K_R(1,1); R_fv=Sol_LS.K_R(2,2); R_s=Sol_LS.K_R(1,2); R_u0=Sol_LS.K_R(1,3); R_v0=Sol_LS.K_R(2,3);
k_c0 = Sol_LS.k_c;
k_d0 = Sol_LS.k_d;

options = optimset('Algorithm',{'levenberg-marquardt', 1e5},'TolX',1e-5,'TolFun',1e-5,  ... 
                     'Jacobian','off','Display','off', 'MaxFunEvals',1300, 'FinDiffType','central');
lb = []; ub = [];

run f_extractParamsForNLS
[X,RESIDUAL] = lsqnonlin(costFunc,X0,lb,ub,options);
run f_composeParamsForNLS

Sol_NLS.R_R_D= R_R_D;
Sol_NLS.R_t_D= R_t_D;
Sol_NLS.K_D= K_D;
Sol_NLS.K_R= K_R;
Sol_NLS.R_P_D = [Sol_NLS.R_R_D*inv(Sol_NLS.K_D), Sol_NLS.R_t_D];
Sol_NLS.k_c = k_c;
Sol_NLS.k_d = k_d;


function F = f_min6pntKinect(X,Sol_LS, f_Opt)

run f_composeParamsForNLS
K_D = K_D;
K_R = K_R;
R_P_D = K_R*[R_R_D R_t_D];
D_P_D = K_D*[eye(3) [0;0;0]];
R_o_e = Sol_LS.R_o_e;
D_o_e = Sol_LS.D_o_e;
Sph = Sol_LS.Sph;

%% Obtain right part of reprojection minimization function either by projecting ellipse center
parfor i=1:length(R_o_e(1,:))
    D_X_vis = [];
    D_X_vis = f_kinect_depth2XYZ(Sph(i).D_U_vis, K_D);
    D_O_s = f_sphereLinLS(D_X_vis);
    beta= 1/(D_O_s(4)^2);
    Q_Ds = [eye(3) - beta*D_O_s([1:3])*D_O_s([1:3])'  -beta*D_O_s([1:3]);
        -beta*D_O_s([1:3])'   -beta];
    D_conic= inv(D_P_D*Q_Ds*D_P_D');
    
    %% Depth;  Fixed Ellipse
    D_o_s_temp = f_ellipse2spherecenter(D_o_e(1:2,i), D_conic, K_D);
    D_o_s(:,i)  = [D_o_s_temp; D_O_s(3)];
    
    %% RGB Ellipse
    R_o_s(:,i) = f_ellipse2spherecenter(R_o_e(:,i), Sph(i).C_R, K_R);
end

XYZ_D = f_depth2XYZ(K_D,D_o_s); %XYZ;%
if f_Opt.case == 3 || f_Opt.case == 5
    XYZ_D(4,:) = [];
    xyz_d = XYZ_D./(ones(3,1)*XYZ_D(3,:));
    [L_rad_d,L_tan_d] = f_calculateDistortion(xyz_d(1:2,:),k_d);
    
    x_d_temp = (xyz_d(1:2,:).*(ones(2,1)*L_rad_d));
    x_d_temp2 = [x_d_temp.*(ones(2,1)*XYZ_D(3,:));
        XYZ_D(3,:)];
    x_temp = [R_R_D, R_t_D]*[x_d_temp2; ones(1,size(x_d_temp2,2))];
else
    x_temp = [R_R_D, R_t_D]*XYZ_D;
end
x_n = x_temp./(ones(3,1)*x_temp(3,:));
% Calculate distortion params:
if f_Opt.case == 3 || f_Opt.case == 5
    [L_rad_r,x_g] = f_calculateDistortion(x_n,k_c);
    % without Tangential Dist:
    x_k = (x_n(1:2,:) .* (ones(2,1)*L_rad_r))+ x_g;
else
    x_k = x_n(1:2,:);
end
% Reconvert in pixels:
u_d = K_R(1:2,1:2)*x_k + K_R(1:2,3)*ones(1,length(x_k));
pix_Err = R_o_s - u_d;

Nmax = max(D_o_s(3,:))+0.25;
Nmin = min(D_o_s(3,:))-0.25;
N = D_o_s(3,:);
S = 1 - exp((Nmax-N)./(Nmin-N));
sumS = sum(S);
W = S/sumS;

F = sqrt(W'.*sum( pix_Err'.^2 ,2)); %W'.*



function F = f_minConicKinect(X,Sol_LS, f_Opt)
run f_composeParamsForNLS
K_D = K_D;
K_R = K_R;
R_P_D = K_R*[R_R_D R_t_D];
D_P_D = K_D*[eye(3) [0;0;0]];
R_o_e = Sol_LS.R_o_e;
D_o_e = Sol_LS.D_o_e;
Sph = Sol_LS.Sph;

%% By sphere fitting (compute sphere center in {D} and projected sphere center in {R}
parfor i=1:length(R_o_e(1,:)),
    D_X_vis = f_kinect_depth2XYZ(Sph(i).D_U_vis, K_D);
    Fit = f_sphereLinLS(D_X_vis);
    beta= 1/(Fit(4)^2);
    Q_Ds = [eye(3) - beta*Fit([1:3])*Fit([1:3])' -beta*Fit([1:3]);
            -beta*Fit([1:3])'   -beta];
    wi=exp(-1/4*(2-D_o_e(3,i))^2);
    left_conic = Sph(i).C_R;
    right_conic= inv(R_P_D*Q_Ds*R_P_D');

    [x0_l, y0_l, a_l, b_l, theta_l, X_l, ax_max, ax_min, resid2] = f_C2param(left_conic);
    [x0_r, y0_r, a_r, b_r, theta_r, X_r, ax_max, ax_min, resid2] = f_C2param(right_conic);
    
    F(i) = wi*norm([x0_l-x0_r; y0_l-y0_r; a_l-a_r; b_l-b_r;(theta_l-theta_r)*pi/180]); %wi*
end
