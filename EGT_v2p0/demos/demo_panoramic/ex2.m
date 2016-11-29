%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 2 - Central Catadioptric Imaging of 3D scene points  
%
close all; clear all; figure(1); hold on; f_3Dwf('k',0.1);

%% Panoramic Camera
H=[rotoy(0)*rotoz(0)*rotox(-180*pi/180) , [-0.2,-0.1,0]';
        0    0    0       ,    1    ];
%quadric=1; a=0.04; b=0.02; r_rim=0.04;%Hyperbola
quadric=2; a=0.03; b=1;    r_rim=0.05; %Parabola
f_3Dpanoramic(H,'g',quadric,a,b,r_rim);     
f_3Dframe(H,'g',0.06,'_{m}');  % Mirror reference frame

%Scene Point
X=[ 0  -.2  .1; 0   .1  .2; .15   .3  .1];
% Camera calibration matrix
K=[10^3   0   320;  0   10^3  240; 0  0  1];
% Point projection
figure(1); axis equal; view(69,18);
[q,Xh]=f_panproj(X,H,K,a,b,quadric,'b:');
title('Example 2 - Central Catadioptric Imaging of 3D scene points');
plot3( X(1,:), X(2,:), X(3,:),'r*'); grid on
f_3Dwfenum(X,'k',0.01);

% Point projection
figure(2); hold on; axis equal; view(0,90);
f_3Dwf('k',0.1); f_3Dpanoramic(H,'g',quadric,a,b,r_rim);     
f_3Dframe(H,'g',0.06,'_{m}'); 
[q,Xh]=f_panproj(X,H,K,a,b,quadric,'b:');
title('Example 2 - Top view for imaging of 3D scene points');
plot3( X(1,:), X(2,:), X(3,:),'r*'); grid on
f_3Dwfenum(X,'k',0.01);

figure(3); plot(q(1,:),q(2,:),'go');
title('Example 2 : CCD Central Catadioptric Camera Plane');
xlabel('u [pixels]'); ylabel('v [pixels]'); grid on; axis equal

