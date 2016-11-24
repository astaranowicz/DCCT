% x1= on the sphere (in <c>)
% xt1= on the virtual plane (in <p>)
% m1 = on the image plane
% R=Rc_m
% t=tc_m (from <matlab> to <camera> frame
function [x1,xt1,m1]=f_catmod(P1,R,t,ell,m,K,col,fig);

tc_m=t; %From o (MATLAB) to C (sphere center)
Rc_m=R;

for i=1:length(P1(1,:)),
%The point P1 in <M> is now expr. in the camera frame <Cc>
    X1([1:3],i) = Rc_m'*P1([1:3],i) - Rc_m'*tc_m;

%Sphere projection (in Cc)
    x1([1:3],i) = X1([1:3],i)/norm(X1([1:3],i)); 
    xs1=x1(1,i);    ys1=x1(2,i);    zs1=x1(3,i); %OK for all

% Projection to the virtual plane
    CcCp=[0; 0; ell]; %As in the fig, Cp is lower than <c>
	xc1([1:3],i)= x1([1:3],i) - CcCp;
	lambda1=(ell-zs1)/(ell+m); % on the triangle equivalence they must be always positive
	xt1([1:3],i)=(1/lambda1)*xc1([1:3],i); %referred to Cp    
    %Projection to the image plane 
    m1([1:3],i)=K*[xt1([1:2],i);1]; %MALIS

end


if fig>=1
% -> PLOT <- %
%%
%% Hereafter all the points must be referred to the MATLAB frame
%% A generic point Pc (belonging to the camera frame <c>) can be mapped to
%% the MATLAB frame <M> as follows:
%%      Pm = Rc_m*Pc + tc_m
%% where Rc_m is "the rotation to move <M> toward <c>
%% and tc_m is the "vector from <M> toward <c>"
%%

%% Parameters - 
% Note that the parameters are changed to their abs values 
% because in the unifying model they must be plotted always positive,
% while in the equivalence with the model they can be negative.
ell=abs(ell);
m=abs(m);
%figure(fig)
set(fig, 'Renderer','OpenGL')
[Xs,Ys,Zs] = sphere(20);
Xs=Xs+tc_m(1);
Ys=Ys+tc_m(2);
Zs=Zs+tc_m(3);
s=surf(Xs,Ys,Zs);
set(s,'FaceAlpha',0.085,'EdgeAlpha',0.095,'FaceColor',col)


%POINTS IN MATLAB
CcM=Rc_m*[0;0;0] + tc_m; %Sphere center
poleM=Rc_m*[0;0;ell] + tc_m; %Sphere pole
PM=Rc_m*[0;0;m]+tc_m; %Plane coordinate
%Frame on the sphere center
Hc=f_Rt2H(Rc_m,tc_m);
Hc=[ Rc_m  tc_m;
    zeros(1,3) 0];
f_3Dframe(Hc,'g',2,'_r');
%Frame of the sphere pole in <M>
Hc=f_Rt2H(Rc_m,poleM); 
Hc=[ Rc_m  poleM;
    zeros(1,3) 0];
f_3Dframe(Hc,'b',.5);
plot3(poleM(1),poleM(2),poleM(3),'k.');
plot3(CcM(1),CcM(2),CcM(3),'k.');

V=[(Rc_m*[1.8;1.8;-abs(m)])'+tc_m';
   (Rc_m*[-1.8;1.8;-abs(m)])'+tc_m';
   (Rc_m*[-1.8;-1.8;-abs(m)])'+tc_m';
   (Rc_m*[1.8;-1.8;-abs(m)])'+tc_m'];

%p=patch(V(:,1),V(:,2),V(:,3),[0 0 0]);
%set(p,'FaceAlpha',0.2,'EdgeAlpha',0.1,'FaceColor',col);
end

if fig==2
for i=1:length(P1(1,:)),
    x1M([1:3],i) =Rc_m*x1([1:3],i) + tc_m; %Point on the sphere
    lambda1b=((ell)-x1(3,i))/(ell+m); % on the triangle equivalence they must be always positive
	xt1b([1:3],i)=(1/lambda1b)*xc1(:,i); %referred to Cp
    xt1M([1:3],i)=Rc_m*(xt1b([1:3],i) + [0;0;ell]) + tc_m; %Point on the virtual plane (must be first expr. in Cc as (xt1+CcCp))
    
    plot3(x1M(1,i),x1M(2,i),x1M(3,i),'g.');
    
         plot3([poleM(1) xt1M(1,i)], [poleM(2) xt1M(2,i)], [poleM(3) xt1M(3,i)],'b:');
         plot3([P1(1,i) CcM(1)],[P1(2,i) CcM(2)],[P1(3,i) CcM(3)],'b:');        
        
    plot3(P1(1,i),P1(2,i),P1(3,i),'r.');
    plot3(xt1M(1,i),xt1M(2,i),xt1M(3,i),'b.');
end

end

end % di tutto il plotta