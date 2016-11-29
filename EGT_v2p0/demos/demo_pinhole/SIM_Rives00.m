%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox  (EGT)  %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Title: "Epipolar Geometry Toolbox application to Visual Servoing"
%  by: Gian Luca Mariottini (Ph.D. Student at SIR Lab- Siena) - January 05
%  (last update)
%
%  From the paper of Patrick Rives "Visual Servoing Based on Epipolar Geometry"
%                    in Proceedings of the 2000 IEEE/RSJ International
%                    Conference on Intelligient Robots and Systems.
%
%  Usage: Simply run this matlab function!
%  -----  
%  Note: If you do not have the Robotics toolbox please set RT (below) = 0
%  ----
%  


%% 1) Initial Setup
clear all
close all
clc

choi=input('Did you have the Robotics Toolbox [by Corke] in your MATLAB? (y/n)>>','s');
if choi=='y',
    RT=1; %if RT=1 ==> Use the Robotics toolbox
else
    RT=0;
end;
figure(1)
hold on
axis equal
grid on

%Scene Points: Random Points  
P=f_3Drandpoint(30,[1 2 .7],1);
  
% %1.1) 3D Scene design: it consists of two rectangular panels
% %Object:panel
% Xpan=1/15*[-5 5 5 -5];
% Ypan=1/15*[ 0 0 0  0];
% Zpan=1/15*[-5 -5 5 5];
% %Orientation first panel
% ang1=0;
% for i=1:length(Xpan),
%     Pan1([1:3],i)=rotoz(ang1)*[Xpan(i);Ypan(i);Zpan(i)];
% end;
% Xp1=Pan1(1,:);
% Yp1=Pan1(2,:);
% Zp1=Pan1(3,:);
% %Orientation second panel
% ang2=pi/4;
% for i=1:length(Xpan),
%     Pan2([1:3],i)=rotoz(ang2)*[Xpan(i);Ypan(i);Zpan(i)];
% end;
% Xp2=Pan2(1,:);
% Yp2=Pan2(2,:);
% Zp2=Pan2(3,:);
% %Orientation third panel
% ang3=-pi/4;
% for i=1:length(Xpan),
%     Pan3([1:3],i)=rotoz(ang3)*[Xpan(i);Ypan(i);Zpan(i)];
% end;
% Xp3=Pan3(1,:);
% Yp3=Pan3(2,:);
% Zp3=Pan3(3,:);
% %Translation of all panels
% n1=length(Xp1);
% P1=[Xp1+0*ones(1,n1) ; Yp1+2*ones(1,n1) ; Zp1-0*ones(1,n1)]; %350 su Y  
% n2=length(Xp2);
% P2=[Xp2-1.1*ones(1,n2) ; Yp2+1.65*ones(1,n2)  ; Zp2+0*ones(1,n2)];   
% n3=length(Xp3);
% P3=[Xp3+0*ones(1,n3) ; Yp3+13*ones(1,n3)  ; Zp3+0*ones(1,n3)];   
% P=[P1 P2];
 for k=1:length(P(1,:)),
          Po([1:4],k)=[P([1:3],k) ; 1];
 end;



  f_3Dwfenum(P);
  f_scenepnt(P,'r*',1);

% 2) Internal Camera parameters
f=1;
u0=0;
v0=0;
Ka=[f 0 u0;
    0 f  v0;
    0 0  1];
Kd=Ka;

% 3) External Camera parameters
  %3.1) Actual Camera (initial position)
    ang_x=20*pi/180;
    ang_y=-10*pi/180;
    ang_z=10*pi/180;
    Ra=rotoz(ang_z)*rotoy(ang_y)*rotox(ang_x);
    ta=[-.2 -1 0.2]';
  %3.2) Desired Camera (target position)
    Rd=eye(3);
    td=[1.1 0.5 0]';

% 4) Plot of cameras
figure(1)
Ha=f_Rt2H(Ra,ta);
f_3Dframe(Ha,'r');
f_3Dcamera(Ha,'r',0.1);
Hd=f_Rt2H(Rd,td);
f_3Dcamera(Hd,'g',0.1);
f_3Dframe(Hd,'r');

% 5-Optional) Definition and plot of Puma 560 6DOF robot arm
if RT==1,
    figure(1);
    Puma560;
    Tpuma_rot=rotz(ang_z)*roty(ang_y)*rotx(ang_x);
    Tpuma_tra=transl(ta);
    Tpuma=Tpuma_tra*Tpuma_rot; %Rotat.&Tranlsation of End-effector   
    Q=ikine560(p560,Tpuma);
    plot(p560,Q);
    axis equal
    view(60,46)
        zoom(2)
end;

%6) Perspective Projection of Scene points
[ud,vd]=f_perspproj(P,Hd,Kd);
[ua,va]=f_perspproj(P,Ha,Ka);
%7) Added noise on corresponding points
sigma=0.00;%
udn=ud+sigma*randn(1,length(ua(1,:)));
vdn=vd+sigma*randn(1,length(ua(1,:)));
uan=ua+sigma*randn(1,length(ua(1,:)));
van=va+sigma*randn(1,length(ua(1,:)));


Ud=[udn;
    vdn];
Ua=[uan;
    van];

%7) Iniatial Epipolar geometry evaluation
figure(2)
axis equal
hold on
title('Actual Image Plane')
plot(uan,van,'r+')

figure(3)
axis equal
hold on
title('Desired Image Plane')
plot(udn,vdn,'gO')
plot(uan,van,'r.')

%Ep.Geom.Estimation/Computation
F=f_Festim(Ua,Ud,4); %FUndamental Matrix Estimation
%[ea,ed,F]=f_epipole(Ha,Hd,Ka,Kd); %Epipole computation

[la,ld]=f_epipline(Ua,Ud,F,1,2,3,'b'); %Epipolar lines computation and plot

% %Starting configuration
theta_x=ta(1);
theta_y=ta(2);
theta_z=ta(3);
omega_x=ang_x;
omega_y=ang_y;
omega_z=ang_z;


col='r.';
beta=0;
if RT==1,
    puma560;
end;
for j=1:200,%number of iterations
  for i=1:length(P(1,:)),    
    %A- Jacobian Matrix Computation  
    x(i)=Ua(1,i);
    y(i)=Ua(2,i);
    Rw2a=Ha([1:3],[1:3])';
    tw2a=-Rw2a*Ha([1:3],4); 
    projpnt([1:3],i)=[Rw2a,tw2a]*Po([1:4],i); %Expression of each feat.point in each camera frame (see EGT manual)    
    Z(i)=projpnt(3,i);
    Lof(:,:,i)=[ -f/Z(i)     0       x(i)/(f*Z(i))     x(i)*y(i)/f^2       -(1+x(i)^2/f^2)      y(i);
                   0      -f/Z(i)    y(i)/(f*Z(i))      1+y(i)^2/f^2      -(x(i)*y(i))/f^2      -x(i)];   
    LofrotT(:,:,i)=[ x(i)*y(i)/f^2       -(1+x(i)^2/f^2)      y(i);
                     1+y(i)^2/f^2      -(x(i)*y(i))/f^2      -x(i)];                               
    %B- Points-Epip.lines distance  
    lin_epip=F*[Ua(:,i) ; 1];
    a(i)=lin_epip(1);
    b(i)=lin_epip(2);    
    LrotT(i,:)=[a(i) b(i)]*LofrotT(:,:,i); %sara' 8x3
    Lep(:,:,i)=[a(i) b(i)]*Lof(:,:,i);
    e(i)=[Ua(:,i) ; 1]'*lin_epip;
  end;
  
    for i=1:length(P(1,:)),
      Lep_tot(i,:)=Lep(:,:,i);
    end;
  
%--> Hybrid Control Law <--   
    Wp=[zeros(3,3); eye(3)];
    e1=pinv(LrotT)*e';
    if norm(e1)<10^-3,
       beta=0.1;
       col='k+';
    end;
    grad_h=[2*(theta_x(j)-td(1)) 2*(theta_y(j)-td(2)) 2*(theta_z(j)-td(3)) 0 0 0]';
    P_K=(eye(6)-Wp*pinv(Wp));
    e_h=Wp*e1+beta*P_K*grad_h; %Control law
    lambda=0.1; %gain factor
    
    U=-lambda*e_h; %control Input to the robot  
    theta_x(j+1)=U(1)+theta_x(j);
    theta_y(j+1)=U(2)+theta_y(j);
    theta_z(j+1)=U(3)+theta_z(j);
    omega_x(j+1)=U(4)+omega_x(j);
    omega_y(j+1)=U(5)+omega_y(j);
    omega_z(j+1)=U(6)+omega_z(j);

        
% EGT --> Update of visual information after camera motion
    Ra=rotoz(omega_z(j+1))*rotoy(omega_y(j+1))*rotox(omega_x(j+1));
    ta=[theta_x(j+1) theta_y(j+1) theta_z(j+1)]';
    Ha=f_Rt2H(Ra,ta);
    [ua,va]=f_perspproj(P,Ha,Ka);
    udn=ud+sigma*randn(1)*ones(1,length(ud));
    vdn=vd+sigma*randn(1)*ones(1,length(vd));
    uan=ua+sigma*randn(1)*ones(1,length(ua));
    van=va+sigma*randn(1)*ones(1,length(va));
    
    Ua=[uan;van];
    Ud=[udn;vdn];    
    F=f_Festim(Ua,Ud,2);
    %[ea,ed,F]=f_epipole(Ha,Hd,Ka,Kd);                   
  
%Plot of results during the servoing.
    figure(4)
    hold on
    title('Error norm - Epiplar Geometry Toolbox')
    xlabel('Iterations')
    plot(j,norm(e),col)
   
    figure(3)
    axis equal
    hold on
    title('Desired Image Plane - Epipolar Geometry Toolbox')
    xlabel('u')
    ylabel('v')
    plot(ua,va,col);

% Camera motion %%-> Note the figure input to functions are REMOVED
    figure(1)    
    hold on
    f_3Dcamera(Ha,'r',0.1);
    text(Ha(1,4)-0.1,Ha(2,4),Ha(3,4)+0.35,'Actual position')
    f_3Dcamera(Hd,'g',0.1);
    text(Hd(1,4)+0.2,Hd(2,4),Hd(3,4),'Desired position')
    f_3Dframe(Ha,'r:');
    f_3Dframe(Hd,'g:');
    f_3Dwfenum(P);
    f_scenepnt(P,'r*',1);
    view(60,46)
    zoom(2)
    axis equal
    grid on
    
    if RT==1,
        Tpuma_rot=rotz(omega_z(j+1))*roty(omega_y(j+1))*rotx(omega_x(j+1));
        Tpuma_tra=transl(ta);
        Tpuma=Tpuma_tra*Tpuma_rot;    
        Q=ikine560(p560,Tpuma);
        plot(p560,Q);
    end;

    axis auto
    
    %Print for GIF animation
    %if mod(j,3)==0,
    %    figura1=strcat('E:\figura1_',num2str(j));
    %    figura3=strcat('E:\figura3_',num2str(j));
    %    figura4=strcat('E:\figura4_',num2str(j));
    %    print( 1, '-djpeg70', figura1);
    %    print( 3, '-djpeg70', figura3);
    %    print( 4, '-djpeg70', figura4);
    %end;    
    
    pause(0.01)
    if j<200
      hold off 
      plot(0,0,'w')
      clear e
    end
end; %\end{FOR}         

%Final Plot
    figure(5)
    hold on
    plot(omega_x*180/pi,'r')
    plot(omega_y*180/pi,'g:')
    plot(omega_z*180/pi,'b--')
    plot(theta_x-td(1),'r-.')
    plot(theta_y-td(2),'g.')
    plot(theta_z-td(3),'b+')
    legend('\omega_x','\omega_y','\omega_z','v_x','v_y','v_z')
                                     