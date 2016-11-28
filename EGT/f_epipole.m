%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %
%%  Epipolar Geometry Toolbox v1.3 (EGT) %
%%                                       %
%%%%%%%  DII - University of Siena  %%%%%%
% 
%  function [ea,ed,F] = f_epipole(Ha,Hd,Ka,Kd)
% 
%  Input:
%  -----
%      Ha    = "Homogenous matrix of actual camera w.r.t. the EGT wrf" 
%      Hd    = "Homogenous matrix of desired camera w.r.t. the EGT wrf" 
%      Ka,Kd = "camera parameters matrix of actual and desired camera
%               respectively"
%
%  Outputs:
%  -------
%      ea,ed = "epipoles in actual and desired camera frame"
%          F = "fundamental matrix"
%
%  Description: 
%  -----------  
%      This function computes the Epipolar geometry between two pin-hole cameras
% 
%  Example:
% 	clear all; close all;figure(2); hold on; figure(3); hold on;
% 	X=[0 , 5]; Y=[10, 4]; Z=[10,-3]; P=[X;Y;Z];     
% 	Rd=eye(3); td=[0,0,0]'; Hd=f_Rt2H(Rd,td);
% 	Ra=rotoy(-pi/6); ta=[-5,-5,0]'; Ha=f_Rt2H(Ra,ta);
% 	Kd=eye(3); Ka=eye(3);
% 	[ud,vd]=f_perspproj(P,Hd,Kd); 
% 	[ua,va]=f_perspproj(P,Ha,Ka);
% 	
% 	[ea,ed,F]=f_epipole(Ha,Hd,Ka,Kd);
% 	figure(2); grid on
% 	title('EGT- Epipolar Geometry - Actual Image plane and epipolar lines')
% 	plot(ea(1),ea(2),'rO'); text(ea(1)+.05,ea(2),'Epipole')
% 	plot(ua,va,'k*'); text(ua+.05,va,'Feature point')
% 	
% 	figure(3); grid on
% 	title('EGT- Epipolar Geometry - Desired Image plane and epipolar lines')
% 	plot(ed(1),ed(2),'gO'); text(ed(1)+.05,ed(2),'Epipole')
% 	plot(ud,vd,'k*'); text(ud+.05,vd,'Feature point')
%
%
% Author:
%     Gian Luca Mariottini 
% Last update:
%     Sept, 18 -2004
function [ea,ed,F]=f_epipole(Ha,Hd,Ka,Kd)
    %It is necessary to premultiplicate 
    Ra=Ha([1:3],[1:3]);
    %ta=-Ra*Ha([1:3],4);%write the translational vector expressed in the camera frame
    Rd=Hd([1:3],[1:3]);
    ta=Ha([1:3],4);%write the translational vector expressed in the camera frame
    td=Hd([1:3],4);%write the translational vector expressed in the camera frame
    
    Trel=Rd'*(ta-td);
    ROTrel=Rd'*Ra;
    E=(f_skew(Trel)*ROTrel);
    F=inv(Kd)'*E*inv(Ka);
    
    ea=null(F);
    ea=ea/ea(3);
    ed=null(F');
    ed=ed/ed(3);
    