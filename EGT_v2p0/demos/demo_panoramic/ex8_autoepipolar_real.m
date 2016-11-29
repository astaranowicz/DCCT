%%% Unifying model for central catadioptric camera
close all; 
%clear all; clc;  

p=0.2; 
fx=17852.7; fy=17912.0;
u0=616.3; v0=628.2;
K=[fx 0 u0;
   0 fy v0;
   0 0  1];

% given the parameters of a real mirror....
a=0.0334;
ell=1;
m=a-ell;

Id=imread('DesiredPan.bmp');
figure(1)
Ud=f_getpoint(Id,3);

Ic_trasl=imread('CurrentPan_onlytransl.bmp');
%Ic_rot=imread('CurrentPan_rotated.bmp');
figure(2)
Uc_trasl=f_getpoint(Ic_trasl,3);

xs_d=f_backcatmod(Ud,m,ell,K); %Backprojection
xs_c=f_backcatmod(Uc_trasl,m,ell,K); %Backprojection

    %-% Autoepipolar Theory on the sphere
     for k=1:length(xs_d(1,:)),
          nvers(1:3,k) = f_normplane(xs_c(1:3,k),xs_d(1:3,k));
          if (norm(nvers(1:3,k))~=0),
               normalvers(:,k)=nvers(:,k);
          end;
     end;
     normal = normalvers;         
     figure(1)
     for h=1:length(normal(1,:)),  
         passo=1/10; ampiezza=500; quadric=2;
         [C(1:3,1:3,h),flag]  = f_panepipconic( normal(1:3,h),K,'b',a,b,quadric,passo,ampiezza );
         axis([0 K(1,3)*2 0 K(2,3)*2]);
     end
     [u,d,v]=svd(normal');

     
     


