%%
% DLT (6point) algorithm that estimates Kd, Dd, R_R_D, R_t_D
%
% Input - R_U - pixel center of the sphere from RGB camera
%         D_m - pixel center of the sphere from Depth Map
%         Kr - RGB camera calibration matrix
%         Dr - RGB camera distortion matrix
%
% Output - M - structure that contains the estimated Kd, Dd, R, and t
%                 Kd -  depth camera calibration matrix
%                 Dd -  depth camera distortion matrix
%                 R - R_R_D - rotation from RGB camera to Depth Camera
%                 t - R_t_D - translation from RGB camera to Depth Camera
%%

function M = f_6pntalgorithmKd(R_U, D_m, Kr, Dr)

XYZ = inv(Kr)*[R_U; ones(1, length(R_U(1,:)))];
XYZ = f_undistort(XYZ, Dr); %TODO: Not necessary at this phase

A=[];
S=[];
%Builds the A matrix used in SVD
for i=1:length(R_U(1,:)),
    D_mZ(:,i) = [ D_m([1:2],i)*D_m(3,i); D_m(3,i); 1];
    A_add = [ zeros(1,4)                -D_mZ(:,i)'           XYZ(2,i)*D_mZ(:,i)' ;
        D_mZ(:,i)'                zeros(1,4)          -XYZ(1,i)*D_mZ(:,i)' ];
    
    A = [  A;
        A_add];
end
[AA,BB,V] = svd( A'*A );
x=V(:,end);
h=x;

R_P_D = [ h(1:4)' ;
    h(5:8)' ;
    h(9:12)'];%OK

[R_R_D, Kd_inv] = qr(R_P_D([1:3],[1:3]));

%% Checks signs on the QR factorization
Kd_inv = Kd_inv/Kd_inv(3,3);
changes=0;
for i=1:3,
    if Kd_inv(i,i)<0;
        Kd_inv(i,:) = - Kd_inv(i,:); % if 1/fu or 1/fv are negative, then change the sign of that row
        R_R_D(:,i) = - R_R_D(:,i); % and change the corresponding column in R
        changes = changes + 1;
    end
end
% Final change of sign, in case only one column of R was changed of sign.
if (mod(changes,2)~=0)||(det(R_R_D) <0),
    R_R_D = -R_R_D;
end

Kd_hat = inv(Kd_inv);
D_t_R = -inv( R_P_D(1:3,1:3)*Kd_hat )*R_P_D(:,4);
R_t_D = -R_R_D*D_t_R;

M.Kd = Kd_hat;
%TODO: Pass initial DCCT_variables.Dd instead
M.Dd = zeros(1, 5);
M.t = R_t_D;
M.R = R_R_D;


end

