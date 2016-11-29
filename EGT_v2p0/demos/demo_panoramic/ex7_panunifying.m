%%% Unifying model for central catadioptric camera
close all; clear all; clc;  
  
P1=[  1, 10,  6 ;
     12,  2, 10 ;
      3,  3,  3];

f=9000; p=0.2; u0=640; v0=480;
% given the parameters of a real mirror....
a=0.03; 
b=1;
% then compute the parameters of the unifying model...

% - PARABOLOID - %
ell=1;
m=a-ell;

K=[f 0 u0;
   0 f v0;
   0 0  1];

theta0=45;
theta_fin=45;
for i=1:theta_fin,
	f1=figure(1); axis equal
	hold on
	%-% Desired Catadioptric Camera
		tc_m1=[0; 0; 0]; %From o (MATLAB) to C (sphere center)
		Rc_m1=eye(3); col1=[.8 .8 .8];
		[xs1,xt1,m1]=f_catmod(P1,Rc_m1,tc_m1,ell,m,K,col1,f1);
    %-% Additive white gaussian noise %-%
        stddev=1;
        for pp=1:length(m1(1,:)),
            m1(:,pp) = m1(:,pp) + [randn(2,1)*stddev/3;0];
        end
        xs1=f_backcatmod(m1,m,ell,K); %Backprojection
	%-% Standard Model - desired frame
		Hpar = [Rc_m1  , tc_m1;
                zeros(1,3) , 1];
		quadric=2; r_rim=0.05; %Parabola
		[qpar1,Xh1]=f_panproj(P1,Hpar,K,a,b,quadric,'b:');

        
	%-% Current Catadioptric Camera
		tc_m2=[10; 10; 0]; %From o (MATLAB) to C (sphere center)
		Rc_m2=rotoz((-theta0+(i-1))*pi/180); col2=[1 0 0];
		[xs2,xt2,m2]=f_catmod(P1,Rc_m2,tc_m2,ell,m,K,col2,f1);
		xs2=f_backcatmod(m2,m,ell,K); %Backprojection
	%-% Standard Model - current frame
		Hpar2 = [Rc_m2  , tc_m2;
                 zeros(1,3) , 1];
		[qpar2,Xh2]=f_panproj(P1,Hpar2,K,a,b,quadric,'b:');
        view(21,28)
        axis equal
        
	f2=figure(2);
	hold on
	for j=1:3,
        plot(m2(1,j),m2(2,j),'rx')
        plot(m1(1,j),m1(2,j),'go')
	end
	grid on, axis equal

    %-% Autoepipolar Theory on the sphere
     for k=1:length(P1(1,:)),
          nvers(1:3,k) = f_normplane(xs1(1:3,k),xs2(1:3,k));
          if (norm(nvers(1:3,k))~=0),
               normalvers(:,k)=nvers(:,k);
          end;
     end;
     normal = normalvers;
     for h=1:length(normal(1,:)),  
         passo=1/10; ampiezza=500; quadric=2;
         [C(1:3,1:3,h),flag]  = f_panepipconic( normal(1:3,h),K,'b',a,b,quadric,passo,ampiezza );
         axis([0 K(1,3)*2 0 K(2,3)*2]);
     end
     [u,d,v]=svd(normal');
     
     %% Bussola traslazione
     im_guess=K*v(:,3);
     im_guess=((im_guess)/im_guess(3));    
     plot([K(1,3) im_guess(1)],[K(2,3) im_guess(2)],'r');
     %%     
     [u1,d1,v1]=svd((normal'*normal));   
     singvd1(i)=d1(3,3)^2/d1(2,2);   
     
     if i~=theta_fin
	pause(0.001);
	figure(1); clf
	figure(2); clf
    end
end


