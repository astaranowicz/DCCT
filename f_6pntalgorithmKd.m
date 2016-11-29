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
% DLT (6point) algorithm that estimates Kd, R_R_D, R_t_D
%
% Input - R_U - pixel center of the sphere from RGB camera
%         D_m - pixel center of the sphere from Depth Map
%         R_K - RGB camera calibration matrix
%
% Output - M - structure that contains the estimated Kd, R, and t
%                 Kd -  depth camera calibration
%                 R - R_R_D - rotation from RGB camera to Depth Camera
%                 t - R_t_D - translation from RGB camera to Depth Camera
%%

function M = f_6pntalgorithmKd(R_U, D_m, R_K)

R_Un = inv(R_K)*[R_U; ones(1, length(R_U(1,:)))];

A=[];
S=[];
%Builds the A matrix used in SVD
for i=1:length(R_U(1,:)),
    D_mZ(:,i) = [ D_m([1:2],i)*D_m(3,i); D_m(3,i); 1];
    A_add = [ zeros(1,4)                -D_mZ(:,i)'           R_Un(2,i)*D_mZ(:,i)' ;
        D_mZ(:,i)'                zeros(1,4)          -R_Un(1,i)*D_mZ(:,i)' ];
    
    A = [  A;
        A_add];
end
[AA,BB,V] = svd( A'*A );
x=V(:,end);
h=x;

R_P_D = [ h(1:4)' ;
    h(5:8)' ;
    h(9:12)'];%OK

[R_R_D, D_Kinv] = qr(R_P_D([1:3],[1:3]));

%% Checks signs on the QR factorization
D_Kinv = D_Kinv/D_Kinv(3,3);
changes=0;
for i=1:3,
    if D_Kinv(i,i)<0;
        D_Kinv(i,:) = - D_Kinv(i,:); % if 1/fu or 1/fv are negative, then change the sign of that row
        R_R_D(:,i) = - R_R_D(:,i); % and change the corresponding column in R
        changes = changes + 1;
    end
end
% Final change of sign, in case only one column of R was changed of sign.
if (mod(changes,2)~=0)||(det(R_R_D) <0),
    R_R_D = -R_R_D;
end

D_K_hat = inv(D_Kinv);
D_t_R = -inv( R_P_D(1:3,1:3)*D_K_hat )*R_P_D(:,4);
R_t_D = -R_R_D*D_t_R;

M.Kd = D_K_hat;
M.t = R_t_D;
M.R = R_R_D;


end

