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
K=[1000  0   640; 
    0   1000  480;
    0    0    1 ];

% Current Panoramic Camera
for i=1:40,
figure(1)
hold on
theta=+pi/4-(i-1)*pi/19/4;

if i>=20, theta=0;end;

t=[0.2,0.2-0*1/200,0]'; %-i/2/100
R=rotoy(0)*rotoz(theta);%*rotox(-180*pi/180));
Hc=[R            , t;
   0    0    0   , 1];
f_3Dpanoramic(Hc,'b',quadric,a,b,r_rim);     
f_3Dframe(Hc,'b',0.06,'_{c}');
    
% Desired Panoramic Camera
Hd=[eye(3) , [0,0,0]';
    0    0    0  , 1];
f_3Dpanoramic(Hd,'g',quadric,a,b,r_rim);     
f_3Dframe(Hd,'g',0.06,'_{d}');

%% SCENE POINTS: Random Points  
%P=f_3Drandpoint(10,[0.4 1 0.4],0.4);
P1=[ 0.4 , 0.2
     0.4 , 0.5
     0.2 , 0.2];
P2=P1 + [0.2, 0.2 ;
          0 , 0   ;
         0.1, 0.1];
P3=P1 + [0.2, 0.2 ;
         -.1, -.1 ;
         0.3, 0.3];
P4=P1 + [ -0.2 , -0.2;
           0.4 , 0.4 ;
          -0.1 ,-0.1];
for ii=1:length(P1(1,:)), 
    plot3(P1(1,ii),P1(2,ii),P1(3,ii),'r*');  %plot in red to identify better
    plot3(P2(1,ii),P2(2,ii),P2(3,ii),'r*');  %plot in red to identify better
    plot3(P3(1,ii),P3(2,ii),P3(3,ii),'r*');  %plot in red to identify better
    plot3(P4(1,ii),P4(2,ii),P4(3,ii),'r*');  %plot in red to identify better
end  
plot3(P1(1,:),P1(2,:),P1(3,:),'r');
plot3(P2(1,:),P2(2,:),P2(3,:),'r');
plot3(P3(1,:),P3(2,:),P3(3,:),'r');
plot3(P4(1,:),P4(2,:),P4(3,:),'r');

%% Baseline: plot %%
figure(1);
[q_bas ,Xh_bas]=f_panproj(t,Hd,K,a,b,quadric,'m:');


%% Feature point: projection %%
[qd1,Xhd1]=f_panproj(P1,Hd,K,a,b,quadric,'g:');
[q1 ,Xhc1]=f_panproj(P1,Hc,K,a,b,quadric,'b:'); %Pair of points

[qd2,Xhd2]=f_panproj(P2,Hd,K,a,b,quadric,'g:');
[q2 ,Xhc2]=f_panproj(P2,Hc,K,a,b,quadric,'b:');

[qd3,Xhd3]=f_panproj(P3,Hd,K,a,b,quadric,'g:');
[q3 ,Xhc3]=f_panproj(P3,Hc,K,a,b,quadric,'b:');

[qd4,Xhd4]=f_panproj(P4,Hd,K,a,b,quadric,'g:');
[q4 ,Xhc4]=f_panproj(P4,Hc,K,a,b,quadric,'b:');

f_3Dwfenum([P1,P2,P3,P4],'k',0.01);  grid on; view(-82,34); axis equal; 
figure(1);
view(36,26)
title('Example 6: Plot of lines into circles in paracatadioptric');

figure(3); hold on;
plot(q1(1,:),q1(2,:),'r*');

figure(3); hold on;
plot(qd1(1,:),qd1(2,:),'gO');
xlabel('u [pixels]'); ylabel('v [pixels]'); grid on;
axis equal

plot(q_bas(1),q_bas(2),'mo'); text(q_bas(1),q_bas(2)+3,'epipole')

f_3d=figure(1);
f_image_d=figure(3);
col='g';
[cd_u1,cd_v1,rd1,nd1,Pd_1]=f_panparconic(P1,Hd,K,a,b,quadric,col,f_3d,f_image_d); 
%P_d1 is in the camera frame, nd_1 is the normal to int. plane 
[cd_u2,cd_v2,rd2,nd2,Pd_2]=f_panparconic(P2,Hd,K,a,b,quadric,col,f_3d,f_image_d);
[cd_u3,cd_v3,rd3,nd3,Pd_3]=f_panparconic(P3,Hd,K,a,b,quadric,col,f_3d,f_image_d);
[cd_u4,cd_v4,rd4,nd4,Pd_4]=f_panparconic(P4,Hd,K,a,b,quadric,col,f_3d,f_image_d);

f_image_c=figure(3);
col='r';
[cc_u1,cc_v1,rc1,nc1,Pc_1]=f_panparconic(P1,Hc,K,a,b,quadric,col,f_3d,f_image_c);
[cc_u2,cc_v2,rc2,nc2,Pc_2]=f_panparconic(P2,Hc,K,a,b,quadric,col,f_3d,f_image_c);
[cc_u3,cc_v3,rc3,nc3,Pc_3]=f_panparconic(P3,Hc,K,a,b,quadric,col,f_3d,f_image_c);
[cc_u4,cc_v4,rc4,nc4,Pc_4]=f_panparconic(P4,Hc,K,a,b,quadric,col,f_3d,f_image_c);


hd1=Pd_1([1:3],1)-Pd_1([1:3],2);
dd1=Pd_2([1:3],1)-Pd_1([1:3],1);
hc1=Pc_1([1:3],1)-Pc_1([1:3],2);
dc1=Pc_2([1:3],1)-Pc_1([1:3],1);

norm_d1=norm(f_skew(Xhd1([1:3],1))*Xhd1([1:3],2));
norm_d2=norm(f_skew(Xhd2([1:3],1))*Xhd2([1:3],2));
norma=norm((f_skew(Pd_2([1:3],1))*Pd_2([1:3],2)));
normb=norm((f_skew(Pd_1([1:3],1))*Pd_1([1:3],2)));
%OK
(f_skew(Pd_2([1:3],1))*Pd_2([1:3],2))/norma - (f_skew(Pd_1([1:3],1)+dd1)*(Pd_1([1:3],2)+dd1))/norma
%OK
nd2*norma - nd1*normb -  (f_skew(hd1)*dd1)

%nd2*norm_d2 - (nd1 + f_skew(hd1)*dd1)*norm_d2 %PROBLEMA: Dovrebbe far zero ma viene diff da zero!
pause


cu_d1=[cd_u1 ; cd_v1];
cu_c1=[cc_u1 ; cc_v1];
cu_d2=[cd_u2 ; cd_v2];
cu_c2=[cc_u2 ; cc_v2];
cu_d3=[cd_u3 ; cd_v3];
cu_c3=[cc_u3 ; cc_v3];
cu_d4=[cd_u4 ; cd_v4];
cu_c4=[cc_u4 ; cc_v4];


e_21 =(cu_c2-cu_c1);
ep_21=(cu_d2-cu_d1);
e_31 =(cu_c3-cu_c1);
ep_31=(cu_d3-cu_d1);
e_32 =(cu_c3-cu_c2);
ep_32=(cu_d3-cu_d2);
e_41 =(cu_c4-cu_c1);
ep_41=(cu_d4-cu_d1);

y1=inv(K)*([e_21;0])/norm(inv(K)*([e_21;0]));
y2=inv(K)*([e_31;0])/norm(inv(K)*([e_31;0]));
y3=inv(K)*([e_32;0])/norm(inv(K)*([e_32;0]));
y4=inv(K)*([e_41;0])/norm(inv(K)*([e_41;0]));

x1=inv(K)*([ep_21;0])/norm(inv(K)*([ep_21;0]));
x2=inv(K)*([ep_31;0])/norm(inv(K)*([ep_31;0]));
x3=inv(K)*([ep_32;0])/norm(inv(K)*([ep_32;0]));
x4=inv(K)*([ep_41;0])/norm(inv(K)*([ep_41;0]));

%% Procrustes Problem Y=RX s.t. R'R=I
B=[x1 x2 x3 x4]';
A=[y1 y2 y3 y4]';
[U,S,V]=svd(A'*B);
R_est=U*V'
display('Press key to continue')
pause


%% Normalized circle centers
% cd_x=-(cd_u-K(1,3))/(2*a/2*K(1,1)); %Qua ci va a/2 perche' dai miei conti 
% cd_y=-(cd_v-K(2,3))/(2*a/2*K(1,1));
% cc_x=-(cc_u-K(1,3))/(2*a/2*K(1,1));
% cc_y=-(cc_v-K(2,3))/(2*a/2*K(1,1));

cd_x1=-(cd_u1-K(1,3))/(2*a/2*K(1,1)); %Qua ci va a/2 perche' dai miei conti 
cd_y1=-(cd_v1-K(2,3))/(2*a/2*K(1,1));
cc_x1=-(cc_u1-K(1,3))/(2*a/2*K(1,1));
cc_y1=-(cc_v1-K(2,3))/(2*a/2*K(1,1));

cd_x2=-(cd_u2-K(1,3))/(2*a/2*K(1,1)); %Qua ci va a/2 perche' dai miei conti 
cd_y2=-(cd_v2-K(2,3))/(2*a/2*K(1,1));
cc_x2=-(cc_u2-K(1,3))/(2*a/2*K(1,1));
cc_y2=-(cc_v2-K(2,3))/(2*a/2*K(1,1));


%% 3-D Rotation Estimation
x1=[ep_ij;0]/norm(ep_ij);
y1=[e_ij;0]/norm(e_ij);

x2=[ep_ik;0]/norm(ep_ik);
y2=[e_ik;0]/norm(e_ik);
% y3=[e_iz;0]/norm(e_iz);
% y4=[e_jk;0]/norm(e_jk);
% y5=[e_jz;0]/norm(e_jz);
% y6=[e_kz;0]/norm(e_kz);


% x3=[ep_iz;0]/norm(ep_iz);
% x4=[ep_jk;0]/norm(ep_jk);
% x5=[ep_jz;0]/norm(ep_jz);
% x6=[ep_kz;0]/norm(ep_kz);
    %%  
    %% R=R{z,phi}*R{y,theta}*R{x,psi}
    ex1=x1(1);
    ex2=x1(2);
    exp1=x2(1);
    exp2=x2(2);
    uai1=y1(1);
    uai1p=y2(1);
    
    Ai=null([-ex1    ex2
            -exp1  exp2]);
    theta_est3D=atan(Ai(1))*180/pi
    psi_est3D=asin(Ai(2))*180/pi
    
    Bi=null([(-cos(theta_est3D)*ex1 - sin(theta_est3D)*sin(psi_est3D)*ex2)  cos(psi_est3D)*ex2 uai1
            (-cos(theta_est3D)*exp1 - sin(theta_est3D)*sin(psi_est3D)*exp2)  cos(psi_est3D)*exp2 uai1p]);
    
        phi_est3D=atan2(Bi(2),Bi(1))*180/pi
        pause
    
    
%%
%% Calibration from Disparity Circles: Image Center
%% (see my papers)
    figure(3);
    dd12=norm(cu_d2-cu_d1);
    ad12=(rd1^2 - rd2^2 + dd12^2)/(2*dd12);
    hd12=sqrt(rd1^2 - ad12^2);
    pd0 = (cu_d2-cu_d1)*ad12/dd12 + cu_d1;
    pd1 = pd0 + [cd_v2-cd_v1; cd_u1-cd_u2]*hd12/dd12;
    pd2 = pd0 - [cd_v2-cd_v1; cd_u1-cd_u2]*hd12/dd12;
    plot([pd1(1),pd2(1)],[pd1(2),pd2(2)],'k')
    plot([pd1(1)],[pd1(2)],'k.')
    plot([pd2(1)],[pd2(2)],'k.')
    
    
    d12=norm(cu_c2-cu_c1);
    a12=(rc1^2 - rc2^2 + d12^2)/(2*d12);
    h12=sqrt(rc1^2 - a12^2);
    p0 = (cu_c2-cu_c1)*a12/d12 + cu_c1;
    p1 = p0 + [cc_v2-cc_v1; cc_u1-cc_u2]*h12/d12;
    p2 = p0 - [cc_v2-cc_v1; cc_u1-cc_u2]*h12/d12;
    plot([p1(1),p2(1)],[p1(2),p2(2)],'k')
    plot([p1(1)],[p1(2)],'k.')
    plot([p2(1)],[p2(2)],'k.')
    
    ld12=f_skew([pd1;1])*[pd2;1];
    l12=f_skew([p1;1])*[p2;1];
    
    centro_est=f_skew(ld12)*l12;
    centro_est=centro_est/centro_est(3)
    plot(centro_est(1),centro_est(2),'k*')
    
    2*a/2*K(1,1)*[cd_x1-cd_x2; cd_y1-cd_y2]*ad12/dd12 - 2*a/2*K(1,1)*[cd_x1;cd_y1];
    lambdas_true=ad12/dd12
    lambdas=lambdas_true;
    lambdas_true-(cd_x1^2+cd_y1^2-cd_x2^2-cd_y2^2)/(2*(cd_x1^2+cd_x2^2+cd_y1^2+cd_y2^2-2*cd_x1*cd_x2-2*cd_y1*cd_y2))-1/2
    
    [cd_x1;cd_y1]*(lambdas-1)/lambdas - [cd_x2;cd_y2];
    
    
    [cd_x1;cd_y1]
    
    %pause
    
%%

%% Probabile misura della distanza trale due viste...distanza angolare tra
%% coppie di features su stessa linea (imm. curr. and deisred)
e11 =q([1:2],1)-cu_c;
e12 =q([1:2],2)-cu_c;
ed11=qd([1:2],1)-cu_d;
ed12=qd([1:2],2)-cu_d;
theta_1=f_atan(e11(2),e11(1))*180/pi - f_atan(e12(2),e12(1))*180/pi;
theta_d1=f_atan(ed11(2),ed11(1))*180/pi - f_atan(ed12(2),ed12(1))*180/pi;
THETA(i)=theta_d1-theta_1;



A=[ep_ij(1) ep_ij(2);
   ep_ij(2) -ep_ij(1)];
y=[e_ij(1);e_ij(2)];
x=inv(A)*y;
x=x/norm(x);
theta_est=atan2(x(2),x(1));
theta;
%figure(3)
%quiver(640,480,[10*cos(theta)],[10*sin(theta)],'b')
%pause
%pause 



%R_R=R*rotox(-180*pi/180);
%[cc_x;cc_y;1]-R*[cd_x;cd_y;1]
%pause

if 1
cd_x1=-(cd_u1-K(1,3))/(2*a/2*K(1,1)); %Qua ci va a/2 perche' dai miei conti 
cd_y1=-(cd_v1-K(2,3))/(2*a/2*K(1,1));
cc_x1=-(cc_u1-K(1,3))/(2*a/2*K(1,1));
cc_y1=-(cc_v1-K(2,3))/(2*a/2*K(1,1));

cd_x2=-(cd_u2-K(1,3))/(2*a/2*K(1,1)); %Qua ci va a/2 perche' dai miei conti 
cd_y2=-(cd_v2-K(2,3))/(2*a/2*K(1,1));
cc_x2=-(cc_u2-K(1,3))/(2*a/2*K(1,1));
cc_y2=-(cc_v2-K(2,3))/(2*a/2*K(1,1));

% cd_x3=-(cd_u3-K(1,3))/(2*a/2*K(1,1)); %Qua ci va a/2 perche' dai miei conti 
% cd_y3=-(cd_v3-K(2,3))/(2*a/2*K(1,1));
% cc_x3=-(cc_u3-K(1,3))/(2*a/2*K(1,1));
% cc_y3=-(cc_v3-K(2,3))/(2*a/2*K(1,1));
end;

%% 3-D rotation estimation
A=[cd_x   0   0    cd_y   0   0    1 0 0;
    0   cd_x  0     0   cd_y  0    0 1 0;
    0     0  cd_x   0     0  cd_y  0 0 1;
   cd_x1   0   0    cd_y1   0   0    1 0 0;
    0   cd_x1  0     0   cd_y1  0    0 1 0;
    0     0  cd_x1   0     0  cd_y1  0 0 1;
   cd_x2   0   0    cd_y2   0   0    1 0 0;
    0   cd_x2  0     0   cd_y2  0    0 1 0;
    0     0  cd_x2   0     0  cd_y2  0 0 1];
%     cd_x3   0   0    cd_y3   0   0    1 0 0;
%     0   cd_x3  0     0   cd_y3  0    0 1 0;
%     0     0  cd_x3   0     0  cd_y3  0 0 1];
y=[cc_x;cc_y;1;cc_x1;cc_y1;1;cc_x2;cc_y2;1];
[U1,S1,V1]=svd(A);
% S_pinv=zeros(6,6);
% for i=1:6
%     S_pinv(i,i)=(S1(i,i))^-1;
% end
% x_s=V1*[S_pinv zeros(6,3);
%         zeros(3,9)]*U1'*y

A_pinv_M=pinv(A,10^-5);
x_s_M=A_pinv_M *y;
pippo=f_unstacked(x_s_M,3);
%% Enforcing the ORTHOGONALITY CONSTRAINT
[U_d,S_d,V_d]=svd(pippo);
R_est_d=U_d*eye(3)*V_d';
%%




%% 2-D rotation estimation
A=[cd_x cd_y;
    cd_y -cd_x];
y=[cc_x;cc_y];
x=pinv(A)*y;
x=x/norm(x);
theta
theta_est=atan2(x(2),x(1))


pause(0.1)
figure(1)
clf
figure(2)
clf
figure(3)
clf
end
clf

