%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox  (EGT) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DemoEGT of Epipolar Geometry Toolbox
% ------------------------------------
%
% In this demo are shown some functions to project points
% on the camera frame <c> after a reflection on two mirrors
%
% Last Update - May 08
%       by Stefano Scheggi, Gian Luca Mariottini
%

    close all; clear all;
    format long e;
    clc;

% Scene Points: Random Points  
    P = f_3Drandpoint(20,[.5 -1 .0],.5);

% Mirror 1 (planar)
    n1 = [0 -1 0]';
    n1 = n1 /norm(n1);
    ds1 = 0;

% Mirror 2 (planar)
    angle = pi/4;
    n2 = [sin(angle) -cos(angle) 0]';
    n2 = n2/norm(n2);
    ds2 = 0;

% Draws plane 
    figure(1); hold on; axis equal; grid on
    plane1 = f_3Dplane(n1,ds1,[0,3],[0,2],[-1.5,1.5]);
    plane2 = f_3Dplane(n2,ds2,[-8,-2],[0,-2],[-1.5,1.5]);

    f_3Dwf('k',1.5); % World reference frame
    f_3Dwfenum(P,'k',0.03);
    f_scenepnt(P,'b*',1); %EGT->Plot of scene pnts.
% Calibration Matrix
    u0 = 640; v0 = 480;
    K = [700 0 u0; 
         0 700 v0; 
         0   0  1];
  
    t = [0.6,-2.8,-1]';
    ang_x = 0;
    ang_y = -23;
    ang_z = 0;
% Rotation Matrix
    R = rotox(ang_x*pi/180)*rotoy(ang_y*pi/180)*rotoz(ang_z*pi/180);
    H = f_Rt2H(R,t);
    f_3Dframe(H,'r',1,'_{d}');
    f_3Dcamera(H,'r',0.3,2); 
    
% Camera reflected by mirror 1
    Hv_m1 = f_reflcamera(n1,ds1,H,0,'g:',0.6,0.15,'_{v_{1}^{[1]}}');
    f_3Dframe(Hv_m1,'g',1,'_{v_{1}^{[1]}}');
    f_3Dcamera(Hv_m1,'g',0.3,2);
% Camera reflected by mirror 2
    Hv_m2 = f_reflcamera(n2,ds2,H,0,'g:',0.6,0.15,'_{v_{1}^{[2]}}');
    f_3Dframe(Hv_m2,'g',1,'_{v_{1}^{[2]}}');
    f_3Dcamera(Hv_m2,'g',0.3,1);
    
% Point Reflected onto mirror 1 
    [reflPoint_1, pntMirr_1, str_pntProj_1, str_pntMirr_1] = f_reflpersproj(n1,ds1,H,K,P); 
    f_scenepnt(pntMirr_1,'g*',1);
    f_3Dwfenum(pntMirr_1,'k',0.03);
    
% Reflected Point onto Mirror 2 
    [reflPoint_2, pntMirr_2, str_pntProj_2, str_pntMirr_2] = f_reflpersproj(n2,ds2,H,K,P); 
    f_scenepnt(pntMirr_2,'r*',1);
    f_3Dwfenum(pntMirr_2,'k',0.03);
    view(-15,66)
    
% Image Planes:plot of superimposed mirrored points
    figure
    hold on
    plot(reflPoint_1(1,:),reflPoint_1(2,:),'r+')
    plot(reflPoint_2(1,:),reflPoint_2(2,:),'gO')


%% EPIPOLAR GEOMETRY COMPUTATION and EPIPOLE COMPARISON
% Ground truth epipolar geom. between <v1> and <v2> 
    [e1_gt,e2_gt,F_gt]=f_epipole(Hv_m1,Hv_m2,K,K);
       
% Epipolar geometry computed from images (in real camera but mirr. points)
  % F matrix estimation (Normalized 8-point algorithm)
      F = f_Festim(reflPoint_1,reflPoint_2,2);
        if F~=zeros(3,3),
           e1=null(F);
           e2=null(F');
           e1=e1/e1(3); e2=e2/e2(3);
        else
           display('EGT warning: EGT can not find epipoles. infty value is substituted');
           e1=[inf inf 1]'; e2=[inf inf 1]';
        end
  
% Epipole comparison: Epipole estimation error
  display(' EGT demo result: Errors in the epipole computation between ground truth and estimated: ')
  epipole_e1_error = norm(e1-e1_gt)
  epipole_e2_error = norm(e2-e2_gt)
  
  
  