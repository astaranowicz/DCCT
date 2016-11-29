%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 3 - Epipolar Geometry Estimation 
%             and EPipolar Conic plot for Central Catadioptric Cameras
%
close all; clear all; 
figure(1); hold on; f_3Dwf('k',0.1);

%quadric=1; a=0.04; b=0.02; r_rim=0.04; % Hyperbola
quadric=2; a=0.03; b=1;    r_rim=0.05; % Parabola

% Camera Calibration Matrix
K=[10^3  0   640; 0  10^3  480;  0    0    1 ];

% Current Panoramic Camera
t=[0.2,0.2,0]'; R=rotoz(pi/4);
Hc=[rotoy(0)*rotoz(pi/4)*rotox(-180*pi/180) , [0.2,0.2,0]';
              0    0    0                   ,       1   ];
f_3Dpanoramic(Hc,'b',quadric,a,b,r_rim);     
f_3Dframe(Hc,'b',0.06,'_{c}');
    
% Desired Panoramic Camera
Hd=[rotoy(0)*rotoz(0)*rotox(-180*pi/180) , [0,0,0]';
              0    0    0                ,   1   ];
f_3Dpanoramic(Hd,'g',quadric,a,b,r_rim);     
f_3Dframe(Hd,'g',0.06,'_{d}');

%% SCENE POINTS: Random Points  
P=f_3Drandpoint(10,[0.4 1 0.4],0.4);

%% IMAGE FORMATION MODEL: point projection and plot
figure(1);
[qd,Xhd]=f_panproj(P,Hd,K,a,b,quadric,'g:');
[q ,Xhc]=f_panproj(P,Hc,K,a,b,quadric,'b:');
f_3Dwfenum(P,'k',0.01);  grid on; view(-82,34); axis equal; 
for i=1:length(P(1,:)), 
    plot3(P(1,i),P(2,i),P(3,i),'r*');  %plot in red to identify better
end  
title('Example 3: Epipolar Geometry for Central Catadioptric Cameras');
figure(2); hold on;
plot(q(1,:),q(2,:),'r*');
title('Example 3: Current camera (blue) CCD camera plane');
xlabel('u [pixels]'); ylabel('v [pixels]'); grid on;
figure(3); hold on;
plot(qd(1,:),qd(2,:),'gO');
title('Example 3: Desired camera (green) CCD camera plane');
xlabel('u [pixels]'); ylabel('v [pixels]'); grid on;

%% EPIPOLAR GEOMETRY %%
% Epipole computation (3D mirror and 2D image plane) and plot (mirror)
    figure(1);
    [ Xcd_e , e_cd] = f_panepipoles(Hd,Hc,a,b,quadric,a,b,quadric,K,K,'g','b');
    Xw_d1=Xcd_e(:,1); Xw_d2=Xcd_e(:,2); Xw_c1=Xcd_e(:,3); Xw_c2=Xcd_e(:,4);
    e_d1=e_cd(:,1) ; e_d2=e_cd(:,2) ; e_c1=e_cd(:,3); e_c2=e_cd(:,4);
    % Epipole plot (image plane)    
    figure(2); 
    plot(e_c1(1),e_c1(2),'b*'); plot(e_c2(1),e_c2(2),'b*');
    text(e_c1(1),e_c1(2),'Epipole'); text(e_c2(1),e_c2(2),'Epipole');
    figure(3);  
    plot(e_d1(1),e_d1(2),'g*'); plot(e_d2(1),e_d2(2),'g*');
    text(e_d1(1),e_d1(2),'Epipole'); text(e_d2(1),e_d2(2),'Epipole');
    % Epipole Estimation: M_ESTIMATOR by Torr
    M=[Xhd([1:3],:); Xhc([1:3],:)];
    errorestimaE=10^-5; iterazioni=100;
    [E,g]=f_panEestim(M,iterazioni,errorestimaE);

%% EPIPOLAR CONICS
passo=1/10; ampiezza=500; figure(2); hold on; npunti=length(Xhd(1,:));
for i=1:npunti,
    figure(2); 
    nvers(:,i) = E'*Xhd([1:3],i); %... in current image ref. frame
    [A2,flag]  = f_panepipconic( nvers(:,i),K,'b',a,b,quadric,passo,ampiezza);   
    figure(3); 
    nversp(:,i)= E*Xhc([1:3],i);  %... in desired image ref. frame
    [A2p,flagp] = f_panepipconic(nversp(:,i),K,'b',a,b,quadric,passo,ampiezza);
end