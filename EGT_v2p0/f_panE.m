%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% f_panE(H1,H2) returns the essential matrix E=R*S for panoramic camera with
% hyperbolic mirror, where R is the rotation matrix between two camera and 
% S is the skew-symmetric matrix of translation vector between two camera
% 
% H1 = Homogeneous matrix containing rotation R respect to EGT frame and 
%      translation vector t respect to world frame of first camera
%     
%
%
% Eleonora Alunno - November 2003

function [E] = f_panE(H1,H2);

if nargin==1
    R12 = H1(1:3,1:3);
    t12 = H1(1:3,4);
elseif nargin==2
    Rwf2mir1 = rotox(-pi);
    Rwf2egt = eye(3);
    Regt2mir1 = H1(1:3,1:3);
    Rwf2mir1 = Rwf2egt*Regt2mir1;
    twf2mir1 = H1(1:3,4);
    Regt2mir2 = H2(1:3,1:3);
    Rwf2mir2 = Rwf2egt*Regt2mir2;
    twf2mir2 = H2(1:3,4);
    R12 = Rwf2mir1'*Rwf2mir2;
    t12 = Rwf2mir1'*(twf2mir2-twf2mir1);
end 
  E = R12*f_skew(t12);
  

%%%%%%rivedere