%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Visual Servoing - Demo #1  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Initial Setup
clear all
close all
clc

%% 1) Scene Points Creation 
    % Object: panel
    Xpan=[-.25  .25 .25 -.25];
    Ypan=[ 0 0 0  0];
    Zpan=[-.25 -.25 .25  .25];
    %Orientation panel
    ang1=0;
    for i=1:length(Xpan),
        Pan1([1:3],i)=rotoz(ang1)*[Xpan(i);Ypan(i);Zpan(i)];
    end;
    Xp1=Pan1(1,:);
    Yp1=Pan1(2,:);
    Zp1=Pan1(3,:);
    n1=length(Xp1);
    P1=[Xp1+0*ones(1,n1) ; Yp1+2*ones(1,n1) ; Zp1-0*ones(1,n1)];
    P=P1;
    n=length(P(1,:));
    
    % Object: random
    %P=f_3Drandpoint(4,[1 2 .7],1);
  
    for k=1:length(P(1,:)),
      Po([1:4],k)=[P([1:3],k) ; 1];
    end;

    f_3Dwfenum(P);
    f_scenepnt(P,'b*',1);

%% 2) Internal Camera parameters
    f=1000;
    u0=1024/2;
    v0=768/2;
    K=[f 0 u0;
        0 f v0;
        0 0  1];
    Kt=K;

%% 3) External Camera parameters
  %3.1) Current Camera (initial position)
    ang_x=-45*pi/180;
    ang_y=0*pi/180;
    ang_z=0*pi/180;
    R = rotoz(ang_z)*rotoy(ang_y)*rotox(ang_x);
    t = [0, 1.3, 1]';
  %3.2) Target Camera (target position)
    R_t = eye(3);
    t_t = [0, 1, 0]';

%% 4) Plot of cameras
    figure(1)
    hold on

    H=f_Rt2H(R,t);
    f_3Dcamera(H,'r',0.1,2);

    Ht=f_Rt2H(R_t,t_t);
    f_3Dcamera(Ht,'g',0.1,2);

RT=1;

% 5-Optional) Definition and plot of Puma 560 6DOF robot arm
if RT==1,
    axis equal
    view(60,46)
    Puma560;
    Tpuma_rot=rotz(ang_z)*roty(ang_y)*rotx(ang_x);
    Tpuma_tra=transl(t);
    Tpuma=Tpuma_tra*Tpuma_rot; %Rotat.&Tranlsation of End-effector   
    Q=ikine560(p560,Tpuma);
    plot(p560,Q);
end;

axis equal
axis auto

%6) Perspective Projection of Scene points
Ut=f_perspproj(P,Ht,Kt);
U  =f_perspproj(P,H,K);

%7) Added noise on corresponding points
sigma=0.00;
Un  = U  + sigma*(randn(2,n));
Utn = Ut + sigma*(randn(2,n));

[ut]=[Utn(1,:)];
[vt]=[Utn(2,:)];
[u]=[Un(1,:)];
[v]=[Un(2,:)];

%7) Iniatial Epipolar geometry evaluation
figure(2)
axis equal
hold on
title('Actual Image Plane','Fontsize',14)
plot(u,v,'r+')
axis([0 1024 0 768])

figure(3)
axis equal
hold on
title('Desired Image Plane','Fontsize',14)
plot(ut,vt,'gO')
axis([0 1024 0 768])

% col='r.';
% beta=0;
% if RT==1,
%     puma560;
% end;

inv_K=inv(K);

i=1;
for j=1:200,%number of iterations
    
    tc_w = t;
    Rc_interm = R;
    Rinterm_w = rotox(-pi/2);
    Rc_w = Rc_interm*Rinterm_w;
    Rw_c =Rc_w';
    
    %-A- Compute the features in the current camera frame and the distance Z(i)
    U_back = inv_K*[U;ones(1,n)];
    x = U_back(1,:);
    y = U_back(2,:);
    Pc = Rc_w'*P-Rc_w'*tc_w*ones(1,n);
    Z = Pc(3,:);
    
    %-B- Interaction Matrix Computation 
    for i=1:n,
      Lx([1:2],[1:6],i) = [ -1/Z(i)   0   x(i)/Z(i) x(i)*y(i)  -(1+x(i)^2) y(i);
                               0     -1/Z  x(i)/Z(i) 1+y(i)^2   -x(i)*y(i) -x(i)];
    end                    
    
    
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
                                     