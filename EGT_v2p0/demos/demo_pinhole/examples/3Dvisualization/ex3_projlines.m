%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 3 - Projection of Lines  for perspective cameras
%
    close all; 
    clear all; 

% Camera calibration parameters
    K=[10^3  0   640;
        0  10^3  480;
        0    0    1 ];

% Figure settings
    figure(1); axis equal; view(43,24); hold on;
% Camera positioning and plot
    Rrpy_a=rotoz(0)*rotoy(0)*rotox(0); to=[-2,0,0]';
    H=f_Rt2H(Rrpy_a,to);    
    f_3Dframe(H,'g:',1); f_3Dcamera(H,'g',.2);        

%-% lINE IN 3D-space %-%
    % First line
    P1=[-2,4,2]';
    P2=[-6,4,0]';
    % Second line
    P4=[2,4,2]';
    P5=[3,4,-2]';
    
    P =[P1 , P2];
    Q= [P4 , P5];
    
    f_scenepnt(P,'r*',1); grid on
    f_3Dwfenum(P);
    plot3([P1(1),P2(1)],[P1(2),P2(2)],[P1(3),P2(3)],'r'); %line passing through that points;         
    ext_a=f_perspline(P,H,K,1,1);
    Xh2a=ext_a(:,3);
    
    figure(2); hold on;
    [ua,va]=f_perspproj(P,H,K);
    plot([ua(1),ua(2)],[va(1),va(2)],'g');
    grid on;