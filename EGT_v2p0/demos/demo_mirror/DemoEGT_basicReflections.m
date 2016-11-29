%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox  (EGT) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DemoEGT of Epipolar Geometry Toolbox
% ------------------------------------
%
% The positions of real and virtual points and real
% and virtual camera, are shown.
%
% Last Update _ March 08
%    by Stefano Scheggi, Gian Luca Mariottini
%

    close all; clear all;
    format long e;
    clc;

% Scene Points: Random Points  
    P = f_3Drandpoint(20,[2.8 -2.0 .0],1);

% Mirror 1  
    n1 = [0 -1 0]';
    n1 = n1 /norm(n1);
    ds1 = 0; % a positive distance "ds1" means oriented AS the normal vector "n1"
    
% Scene Points: reflected points (and plot)
    figure(1); hold on; axis equal; grid on
    P_rifl_mirr = f_reflpoint(n1,ds1,P,1,'g*');
    f_3Dwfenum(P_rifl_mirr,'k',0.03);
    
% Draw plane 
    plane1 = f_3Dplane(n1,ds1,[-2,4],[0,2],[-3,1]);
    f_3Dwf('k',1.5); % World reference frame
    f_3Dwfenum(P,'k',0.03);
    f_scenepnt(P,'r*',1); %EGT->Plot of scene pnts.
    
% Calibration Matrix
    u0 = 640; 
    v0 = 480;
    K = [700   0 u0; 
           0 700 v0; 
           0   0  1];

% Camera <c> pose      
    t = [0.6,-2.8,-.5]';% Camera translation (from the <m> frame to the camera <c>)  
    ang_x =  0;
    ang_y = 23;
    ang_z =  0;
    R = rotox(ang_x*pi/180)*rotoy(ang_y*pi/180)*rotoz(ang_z*pi/180);    % Camera Rotation (see EGT manual for convention)   
    H = f_Rt2H(R,t);
    % ... and camera visualization 
   
    f_3Dframe(H,'r',1,'_{c}');
    f_3Dcamera(H,'r',0.3,2); 

% Virtual Camera (Reflected Camera)
    K_v=K;
    H_v = f_reflcamera(n1,ds1,H,0,'g:',0.6,0.15,'_{v_{1}^{[1]}}');
    f_3Dcamera(H_v,'g',0.3,2);
    f_3Dframe(H_v,'g',1,'_{v_{1}^{[1]}}');
    
    view(50,46)
    title(' EGT demo - Mirror reflection of points and camera ');

  