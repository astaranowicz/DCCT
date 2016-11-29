% P=point in the camera frame
% n=normal to the interpretation plane
function [cx,cy,r,n,P]=f_panparconic(P,H,K,a,b,quadric,col,f_3d,f_image);

[q,Xh]=f_panproj(P,H,K,a,b,quadric);
n=(f_skew(Xh([1:3],2))*Xh([1:3],1))/(norm(f_skew(Xh([1:3],2))*Xh([1:3],1))); % Normal vector of the interpretation plane (desired image plane)
P([1:3],1)=H([1:3],[1:3])'*P([1:3],1)-H([1:3],[1:3])'*H([1:3],4); %Feature point referred to the <c> frame
P([1:3],2)=H([1:3],[1:3])'*P([1:3],2)-H([1:3],[1:3])'*H([1:3],4);

c_n=H([1:3],[1:3])*[mean([P(1,:),0]) ; mean([P(2,:),0]); mean([P(3,:),0])]+H([1:3],4); %centroid of the interpretation plane (then referred to the matlab frame, for plotting)
n_M=H([1:3],[1:3])*n;
if (exist('f_3d')~=0)&&(f_3d~=0),
    figure(f_3d);
    hold on;
    quiver3(c_n(1),c_n(2),c_n(3),n_M(1),n_M(2),n_M(3),0.5,col);
end

alpha=K(1,1);
u0=K(1,3);
v0=K(2,3);
c=[u0;v0]-2*(a/2)*alpha*[n(1)/n(3); 
                         n(2)/n(3)];
cx=c(1);
cy=c(2);
r=sqrt(4*(a/2)^2*alpha^2*(1+(n(1)/n(3))^2+(n(2)/n(3))^2));
if (exist('f_image')~=0)&(f_image~=0),
    figure(f_image);
    hold on;
    i=1;
    for T=0:2*pi/360:2*pi,
      x(i)=(r)*cos(T)+cx;
      y(i)=(r)*sin(T)+cy;
      i=i+1;
    end
    plot(x,y,strcat(col,':'));
    plot(cx,cy,strcat(col,'x'));
end