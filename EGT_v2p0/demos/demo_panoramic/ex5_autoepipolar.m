%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 5 - Auto-epipolar Property for Central Catadioptric Cameras
%
% References:  G.L.Mariottini, E.Alunno, J.Piazzi, D.Prattichizzo 
%              

clear all
close all;

figure
text(0.1,.95,'ENLARGE THE FIGURE')
text(0.1,.9 ,'The auto-epipolar property for central catadioptric cameras is based on the construction of "bi-osculating conics".')
text(0.1,.85,'Bi-conics have been introduced by Mariottini,Alunno,Piazzi,Prattichizzo in "Epipole-Based Visual Servoing with Central');
text(0.1,.80,'Catadioptric Camera", submitted at ICRA05.')
text(0.1,.75,'When they all intersect at only one point then it can be show that both cameras have the same orientation.')
text(0.1,.70,'A functional cost based on SVD will be zero in this case (Auto-Epipolar Property).')
text(0.1,.60,'Bi-conics and Auto-Epipolar Property have been used for Visual servoing purposes at the University of Siena')
text(0.1,.50,'Press any KEY to see the animation!')
pause

close all, clear all,


format long
quadric=2; a=0.03; b=1; r_rim=0.05;   %Parabola
%quadric=1; a=0.02; b=0.02; r_rim=0.04;%Hyperbola
      
fig1 = 1;
fig2 = 2;
fig3 = 3;
fig4 = 4;

%% CUBE %%
  %X=10^-1*[-5,  5,  5, -5, -5,  5,  5, -5];
  %Y=10^-1*[-5, -5, -5, -5,  5,  5,  5, 5];
  %Z=10^-1*[1, 1,  5,  5, 1, 1,  5, 5];
  %P=[X;Y;Z];
%% Scene Points: Random Points  
  %P=f_3Drandpoint(4,[0 0.5 1],.7);  
%% Planar %%
  X=10^-1*[-5,  5,  5];%, 1];
  Y=10^-1*[-5, -5, -5];%,-1];
  Z=10^-1*[1,   1,  5];%, 1/5];
  P=[X;Y;Z];
  
  
%% DESIRED CAMERA %%  
Rd = rotoz(0)*rotoy(pi)*rotox(0);
td = [0,0,0]';
Hd = [      Rd       , td ;
         0   0   0   ,  1];
%Scene Points
K=[10^3  0   640;
    0  10^3  480;
    0    0    1 ];

%% NOISE %%
sigma=0;
    
figure(fig3), 
hold on, view(30,42), grid on, axis square
f_3Dpanoramic(Hd,'g',quadric,a,b,r_rim);     
f_3Dframe(Hd,'g',0.06,'_{m}');

[q1,Xhmir1]=f_panproj(P,Hd,K,a,b,quadric,'g:'); %% projection and back-projection (desired camera)
sigma=0;
for i=1:length(q1(1,:)),
     q1(1,i)=q1(1,i)+randn(1)*sigma;
     q1(2,i)=q1(2,i)+randn(1)*sigma;         
end;


Xhmird   = f_daqaXhmir(q1,K,a,b,quadric); %In the mirror frame
Xhmird_m = Rd*Xhmird([1:3],:); %In the MATLAB frame
plot3(Xhmird_m(1,:),Xhmird_m(2,:),Xhmird_m(3,:),'m.');

%% ACTUAL CAMERA %%
ang_in = -50*pi/180;
ang_fin = 50*pi/180;
stepang = 1;
steps = (ang_fin-ang_in)*(180/pi)/stepang;
ta = [-0.2,-0.1,0]'; 

for i=1:steps+1,     
     ang(i) = (ang_in+(i-1)*stepang*pi/180);      
     Ra = rotoy(pi)*rotoz(ang(i)); % "roll-pitch-yaw"
     Ha = [       Ra      , ta;
               0  0  0    ,  1];
     E=f_skew(ta)*Ra;
     fro_norm(i)=norm(E+E','fro');
    
     figure(fig3);
     set(fig3,'Position',[10 200 500 500]);
     title('Example 5 (EGT) - The actual camera rotates until all bi-conics intersect in 1 point');
     hold on;
     view(30,42), grid on;
     f_3Dpanoramic(Hd,'g',quadric,a,b,r_rim);     
     f_3Dframe(Hd,'g',0.06,'_{m}');
     f_3Dpanoramic(Ha,'b',quadric,a,b,r_rim);     
     f_3Dframe(Ha,'b',0.06,'_{m}');
    
     %% projection and back_projection (actual camera)
     [q2,Xhmir2]=f_panproj(P,Ha,K,a,b,quadric,'b:');     
     for c1=1:length(q2(1,:)),
         q2(1,c1)=q2(1,c1)+randn(1)*sigma;
         q2(2,c1)=q2(2,c1)+randn(1)*sigma;         
     end;     

     Xhmira = f_daqaXhmir(q2,K,a,b,quadric);
     Xhmira_m=Ha*Xhmira([1:4],:);
     plot3(Xhmira_m(1,:),Xhmira_m(2,:),Xhmira_m(3,:),'r.');
     
%      % Backprojection on the sphere of radius rho
%      rho=1;
%      for st=1:length(q2(1,:)),
%         f=10^-2;
%         alphax=K(1,1);
%         u0=K(1,3);
%         v0=K(2,3);
%         ua=q2(1,st);
%         va=q2(2,st);
%         ud=q1(1,st);
%         vd=q1(2,st);  
%         betaa=((ua-u0)^2+(va-u0)^2)/(alphax*f);
%         betad=((ud-u0)^2+(vd-u0)^2)/(alphax*f);
%         Xs_a(:,st)=[rho*(ua-u0)/(alphax)*(2/(1+betaa));
%                     rho*(va-v0)/(alphax)*(2/(1+betaa));
%                     rho*(1-betaa)/(1+betaa)];
%         Xs_d(:,st)=[rho*(ud-u0)/(alphax)*(2/(1+betad));
%                     rho*(vd-v0)/(alphax)*(2/(1+betad));
%                     rho*(1-betad)/(1+betad)];
%       end;

% Orthographic Backprojection on the sphere of radius rho
%      rho=10;
%      for st=1:length(q2(1,:)),
%         alphax=K(1,1);
%         u0=K(1,3);
%         v0=K(2,3);
%         ua=q2(1,st);
%         va=q2(2,st);
%         ud=q1(1,st);
%         vd=q1(2,st);  
%         Xs_a(:,st)=[(ua-u0)/(alphax);
%                     (va-v0)/(alphax);
%                     sqrt(rho^2 - ((ua-u0)^2+(va-v0)^2)/(alphax^2)) ];
%         Xs_d(:,st)=[(ud-u0)/(alphax);
%                     (vd-v0)/(alphax);
%                     sqrt(rho^2 - ((ud-u0)^2+(vd-v0)^2)/(alphax^2)) ];
%       end;
     
     
     for k=1:length(q1(1,:)),
          nvers(1:3,k) = f_normplane(Xhmird(1:3,k),Xhmira(1:3,k));
          if (norm(nvers(1:3,k))~=0),
               normalvers(:,k)=nvers(:,k);
          end;
     end;
     normal = normalvers;
    
    passo=1/10; ampiezza=500; 
    for h=1:length(normal(1,:)),
            figure(fig2);
            hold on;
            plot(q1(1,:),q1(2,:),'go');
            plot(q2(1,:),q2(2,:),'r+');       
            legend('Features in current view (rotating camera)','Features in target view (fixed)')   
            [C(1:3,1:3,h),flag]  = f_panepipconic( normal(1:3,h),K,'b',a,b,quadric,passo,ampiezza);
    end

     figure(fig2); 
     set(fig2,'Position',[460 100 600 600]);
     title('Example 5 (EGT) - The bi-conics in the image plane (blue) change as the camera rotates');
     axis equal;
     axis([500 750 300 600]); 
     grid on;     
     
    
    figure(3)
    axis equal
     
    pause(0.01), 
     if ang(i)==0,
            figure(2)
            [Xcd_e,e_cd]= f_panepipoles(Hd,Ha,a,b,2,a,b,2,K,K,'g','r');
            e_d1=e_cd(:,1) ; e_d2=e_cd(:,2) ; e_c1=e_cd(:,3); e_c2=e_cd(:,4);
            % Epipole plot (image plane)    
            figure(2); 
            plot(e_c1(1),e_c1(2),'b*'); plot(e_c2(1),e_c2(2),'b*');
            text(e_c1(1),e_c1(2),'Epipole'); text(e_c2(1),e_c2(2),'Epipole');
            figure(3);  
            plot(e_d1(1),e_d1(2),'g*'); plot(e_d2(1),e_d2(2),'g*');
            text(e_d1(1),e_d1(2),'Epipole'); text(e_d2(1),e_d2(2),'Epipole');
           title('Example 5- Bi-conics intersect at the epipoles when cameras have the same orient.!!!')
           pause(4);         
     end
     if i<steps,
           figure(2),
           hold off, 
           plot(0,0,'w');    
           figure(3),
           hold off, 
           plot(0,0,'w');    
     end;
     r(i)=rank(normal,10^-9);
     [u,d,v]=svd(normal);
     [u1,d1,v1]=svd(normal'*normal);
     singvd(i)=min([d(1,1),d(2,2),d(3,3)]);     
     singvd1(i)=min([d1(1,1),d1(2,2),d1(3,3)]);   
end

figure, hold on
plot(ang*180/pi,(singvd./(norm(singvd))),'b')
plot(ang*180/pi,(singvd1./(norm(singvd1))),'r')
xlabel('Relative orientation between cameras : \theta (degrees)');
ylabel('sv')
title('Singular value behavior - It is zero when cameras have the same orientation')

figure
plot(r)
title('Rank of the bi-conic matrix');