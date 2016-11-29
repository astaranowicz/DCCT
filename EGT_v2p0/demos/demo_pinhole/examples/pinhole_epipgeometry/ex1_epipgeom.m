%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% by Gian Luca Mariottini
%
% Computation of epipolar geometry entities (lines, epipole) and
% visualization.
%
clear all
close all

X=[0 , 5]; Y=[10, 4]; Z=[10,-3]; P=[X;Y;Z];     
figure(2); hold on; figure(3); hold on;
Rd=eye(3); td=[0,0,0]'; Hd=f_Rt2H(Rd,td);
Ra=rotoy(-pi/6); ta=[-5,-5,0]'; Ha=f_Rt2H(Ra,ta);

Kd=eye(3); Ka=eye(3);
Ud=f_perspproj(P,Hd,Kd); 
Ua=f_perspproj(P,Ha,Ka);
ua=Ua(1,:); va=Ua(2,:);
ud=Ud(1,:); vd=Ud(2,:);

[ea,ed,F]=f_epipole(Ha,Hd,Ka,Kd);
figure(2)
title('EGT- Epipolar Geometry - Actual Image plane and epipolar lines')
plot(ea(1),ea(2),'rO'); text(ea(1)+.05,ea(2),'Epipole')
plot(ua,va,'k*'); text(ua+.05,va,'Feature point')
grid on

figure(3)
title('EGT- Epipolar Geometry - Desired Image plane and epipolar lines')
plot(ed(1),ed(2),'gO'); text(ed(1)+.05,ed(2),'Epipole')
plot(ud,vd,'k*'); text(ud+.05,vd,'Feature point')
grid on

[la,ld]=f_epipline(Ua,Ud,F);

