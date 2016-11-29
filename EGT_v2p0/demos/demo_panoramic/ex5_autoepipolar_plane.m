%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 5 - Auto-epipolar Property for Central Catadioptric Cameras
%
% References:  G.L.Mariottini, E.Alunno, J.Piazzi, D.Prattichizzo 
%              
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
a = 3; 
b = 1;
r_rim = 2;
      
fig1 = 1;
fig2 = 2;
fig3 = 3;
fig4 = 4;
 
 %1.1) 3D Scene design: it consists of two rectangular panels
%Object:panel
Xpan=1*[ 0 0   0  0];
Ypan=1*[ 5 5  -5 -5];
Zpan=1*[ 5 -5 -5  5];
%Orientation first panel
ang1=0;
for i=1:length(Xpan),
    Pan1([1:3],i)=rotoy(ang1)*[Xpan(i);Ypan(i);Zpan(i)];
end;
Xp1=Pan1(1,:);
Yp1=Pan1(2,:);
Zp1=Pan1(3,:);
% %Orientation second panel
ang2=0;
for i=1:length(Xpan),
    Pan2([1:3],i)=rotoy(ang2)*[Xpan(i);Ypan(i);Zpan(i)];
end;
Xp2=Pan2(1,:);
Yp2=Pan2(2,:);
Zp2=Pan2(3,:);
%Orientation third panel
ang3=-pi/4;
for i=1:length(Xpan),
    Pan3([1:3],i)=rotoz(ang3)*[Xpan(i);Ypan(i);Zpan(i)];
end;
Xp3=Pan3(1,:);
Yp3=Pan3(2,:);
Zp3=Pan3(3,:);
%Translation of all panels
n1=length(Xp1);
P1=[Xp1-5*ones(1,n1) ; Yp1+0*ones(1,n1) ; Zp1+6*ones(1,n1)]; %350 su Y  
n2=length(Xp2);
P2=[Xp2-10*ones(1,n2) ; Yp2+0*ones(1,n2)  ; Zp2+6*ones(1,n2)];   
n3=length(Xp3);
P3=[Xp3+0*ones(1,n3) ; Yp3+13*ones(1,n3)  ; Zp3+0*ones(1,n3)];   
%Final set of points
P=[P1 P2];


  
%% DESIRED CAMERA %%  
Rd = rotoz(0)*rotoy(0)*rotox(0);
td = [0,0,0]';
Hd = [ rotoy(pi)*Rd  , td ;
         0   0   0   ,  1];
K = eye(3);
    K=[10^1  0   320;
        0  10^1  240;
        0    0    1 ];
%% projection and back-projection (desired camera)
[q1,Xhmir1] = f_panproj(P,Hd,K,a,b);
Xhmird = f_daqaXhmir(q1,K,a,b);

%figure(fig1), hold on, axis equal

figure(fig3), 
hold on, view(30,42), grid on, axis square
f_3Dpanoramic(Hd,'g',1,a,b,r_rim);
f_3Dframe(Hd,'g',3,'_d');



%% ACTUAL CAMERA %%
ang_in = -100*pi/180;
ang_fin = 100*pi/180;
stepang = 5;
steps = (ang_fin-ang_in)*(180/pi)/stepang;

ta = [26, 15, 0]'; 

for i=1:steps+1,     
     ang(i) = (ang_in+(i-1)*stepang*pi/180);      
     Ra = rotoz(ang(i)); % "roll-pitch-yaw"
     Ha = [rotoy(pi)*Ra , ta;
            0  0   0    ,  1];
        
     E=f_skew(ta)*Ra;
     fro_norm(i)=norm(E+E','fro');
      
     %% projection and back_projection (actual camera)
     [q2,Xhmir2] = f_panproj(P,Ha,K,a,b);
     Xhmira = f_daqaXhmir(q2,K,a,b);
      
     for k=1:length(q1(1,:)),
          nvers(1:3,k) = f_normplane(Xhmird(1:3,k),Xhmira(1:3,k));
          if norm(nvers(1:3,k))~=0
               normalvers(:,k)=nvers(:,k);
          end
     end
     normal = normalvers;
    
     
     for h=1:length(normal(1,:)),
           figure(fig2);
           hold on;
           plot(q1(1,:),q1(2,:),'go');
           plot(q2(1,:),q2(2,:),'r+');       
           legend('Features in current view (rotating camera)','Features in target view (fixed)')
           [C(1:3,1:3,h),flag] = f_projectconic(normal(1:3,h),K,fig2,'b',a,b);                      
           if flag~=0,
               disp('  This plane, is the one passing through the mirror reprojection of the image points')
               disp(q1(:,h));
               disp(q2(:,h));     
           end   
     end

if 1,     
%% PLOT {begin}
     figure(fig2), 
     set(fig2,'Position',[460 100 600 600])
     title('Example 5 (EGT) - The bi-conics in the image plane (blue) change as the camera rotates')
     axis([320-2 320+2 240-2 240+2]), 
     grid on;     
     
     figure(fig3),
     set(fig3,'Position',[10 200 500 500])
     title('Example 5 (EGT) - The actual camera rotates until all bi-conics intersect in 1 point')
     hold on
     view(30,42), grid on, 
     f_3Dpanoramic(Hd,'g',1,a,b,r_rim);
     f_3Dframe(Hd,'g',3,'_d');
     f_3Dpanoramic(Ha,'r',1,a,b,r_rim);
     f_3Dframe(Ha,'r',3,'_a');
     f_3Dwfenum(P);  
     plot3(P(1,:),P(2,:),P(3,:),'r*');grid on;
     axis equal 
     
     pause(0.0001), 
     
     if ang(i)==0,
           figure(2)
           [e1d,e2d,e1a,e2a] = f_panepipoles(Hd,Ha,3,1,3,1,K,K,'g.','r.',-1);
           plot(e1a(1),e1a(2),'cs');
           plot(e2a(1),e2a(2),'cs');
           text(e1a(1),e1a(2),'  epipole');
           text(e2a(1),e2a(2),'  epipole'); 
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
%% PLOT {end}
end
     
     
     [u,d,v]=svd(normal);
     singvd(i)=min([d(1,1),d(2,2),d(3,3)]);     
     ugual(i)=d(1,1)-d(2,2);
end

figure, hold on
plot(ang*180/pi,singvd.^2)
plot(ang*180/pi,ugual.^2)
xlabel('Relative orientation between cameras : \theta (degrees)');
ylabel('sv')
title('Singular value behavior - It is zero when cameras have the same orientation')