%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.1 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% [x,y,z,nx,ny,nz]=f_3Dsurface; 
% [x,y,z,nx,ny,nz]=f_3Dsurface(ogg); 
% [x,y,z,nx,ny,nz]=f_3Dsurface(ogg,R); 
% [x,y,z,nx,ny,nz]=f_3Dsurface(ogg,R,t);
% [x,y,z,nx,ny,nz]=f_3Dsurface(ogg,R,t,param);
% [x,y,z,nx,ny,nz]=f_3Dsurface(ogg,R,t,param,n);
% [x,y,z,nx,ny,nz]=f_3Dsurface(ogg,R,t,param,n,disegna);
%
%Syntax:
%-------
%     ogg= "integer representing the type of represented surface.
%           ogg=2 --> "sphere"
%           ogg=1 --> "torus"
%           ogg=3 --> "ellipse"
%           ogg=4 --> "infinity"
%           ogg=5 --> "spinner"
%           ogg=6 --> "cone"
%     t= "traslational vector wrt the wrf."
%     R= "rotation matrix rpy wrt the wrf."
%     param= "geoemtrical parameters of selected surface (variable wrt ogg)"
%                if ogg==2 ==> param=[par1 par2] where r1=internal radius
%                                                      r2=external radius
%                if ogg==1 ==> param="sphere radius"
%
%Example 1:
%       f_3Dsurface(2,rotoy(pi/6),[10,10,0]',[3 1 1],30,1); %torus
%
%Example 2:
%       f_3Dsurface(1,eye(3),[10,10,0]',[3 1],13,1);
%
% Author:
%   Gian Luca Mariottini
% Last Update:
%   Sept, 18 - 2004
function [x,y,z,nx,ny,nz]=f_3Dsurface(ogg,R,t,param,n,disegna)

%SPHERE
if ogg == 1,
    % -pi/2 < theta < pi/2 
    % 0 <= phi < 2*pi 
    passo=param(2);
    rad=param(1);
    theta = ((0:passo:n-1)*pi/(n-1)-pi/2); %in radians
    phi = ((0:passo:n-1)*2*pi/(n-1)); 
    %
    %r(theta,phi)=( rad*cos(theta)*cos(phi) , rad*cos(theta)*sin(phi) , rad*(sin(phi)))
    %
    x_0=rad*cos(theta).'*cos(phi);
    y_0=rad*cos(theta).'*sin(phi);
    z_0=rad*sin(theta).'*ones(1,n);
%     nx=cos(theta).'*cos(phi);
%     ny=cos(theta).'*sin(phi);
%     nz=sin(theta).'*ones(1,n);

end;

%Torus
if ogg == 2,
    %
    % par1=radius of internal circle
    % par2=center of internal circle= (radius -1) of external circle
    % par3=angle_step
    % 
    %r(phi,theta)=( rad*cos(phi) , par2*(sin(phi)+par1)cos(theta) , par2*(sin(phi)+par1)sin(theta))
    %
    passo=param(3);
    phi = ((0:passo:n-1)*2*pi/(n-1)); 
    theta = ((0:passo:n-1)*2*pi/(n-1));
    x_0=param(2)*cos(phi).'*ones(1,n);
    y_0=param(2)*(sin(phi)+param(1)).'*cos(theta);
    z_0=param(2)*(sin(phi)+param(1)).'*sin(theta);
end;

%Ellipse
if ogg == 3,
    %
    % par1 = a (semiaxis alog x_axis)
    % par2 = b (semiaxis alog y_axis)
    % par3 = angle_step
    % 
    %r(phi,theta)=( par1*cos(phi) , par2*sin(phi)*cos(theta) , par2*sin(phi)*sin(theta) )
    %
    passo=param(3);
    phi = ((0:passo:n-1)*2*pi/(n-1)); 
    theta = ((0:passo:n-1)*2*pi/(n-1));
    x_0=param(1)*cos(phi).'*ones(1,n);
    y_0=param(2)*sin(phi).'*cos(theta);
    z_0=param(2)*sin(phi).'*sin(theta);
end;


%Infinity
if ogg == 4,
    %
    % par1 = a (semiaxis alog x_axis)
    % par2 = b (semiaxis alog y_axis)
    % par3 = angle_step
    % 
    %r(phi,theta)=( par1*cos(phi) , par2*sin(phi)*cos(theta) , par2*sin(phi)*sin(theta) )
    %
    passo=param(3);
    phi = ((0:passo:n-1)*2*pi/(n-1)); 
    theta = ((0:passo:n-1)*2*pi/(n-1));
    x_0=param(1)*cos(phi).'*ones(1,n);
    y_0=param(2)*(sin(phi+theta)).'*cos(theta);
    z_0=param(2)*sin(phi+theta).'*sin(theta);
end;

%Spinner %trottola
if ogg == 5,
    %
    % par1 = a (semiaxis alog x_axis)
    % par2 = b (semiaxis alog y_axis)
    % par3 = angle_step
    % 
    %r(phi,theta)=( par1*cos(phi) , par2*sin(phi)*cos(theta) , par2*sin(phi)*sin(theta) )
    %
    passo=param(3);
    phi = ((0:passo:n-1)*2*pi/(n-1)); 
    theta = ((0:passo:n-1)*2*pi/(n-1));
    x_0=param(1)*atan(phi).'*ones(1,n);
    y_0=param(2)*sin(phi).'*cos(theta);
    z_0=param(2)*sin(phi).'*sin(theta);
end;

%Cone
if ogg == 6,
    %
    % par1 = a (semiaxis alog x_axis)
    % par2 = b (semiaxis alog y_axis)
    % par3 = angle_step
    % 
    %r(phi,theta)=( par1*cos(phi) , par2*sin(phi)*cos(theta) , par2*sin(phi)*sin(theta) )
    %
    passo=param(3);
    phi = ((0:passo:n-1)*2*pi/(n-1)); 
    theta = ((0:passo:n-1)*2*pi/(n-1));
    x_0=(param(1))*sinh(phi).'*ones(1,n);
    y_0=(param(1))*cosh(phi).'*cos(theta);
    z_0=(param(1))*cosh(phi).'*sin(theta);
end;

%Rotation & traslation
for i=1:n,
     for j=1:n,
         x_t(i,j) = R(1,1)*x_0(i,j)+R(1,2)*y_0(i,j)+R(1,3)*z_0(i,j)+t(1);
         y_t(i,j) = R(2,1)*x_0(i,j)+R(2,2)*y_0(i,j)+R(2,3)*z_0(i,j)+t(2);
         z_t(i,j) = R(3,1)*x_0(i,j)+R(3,2)*y_0(i,j)+R(3,3)*z_0(i,j)+t(3);
     end;
end;

x=x_t;
y=y_t;
z=z_t;

%Normals Matrix
[nx,ny,nz]=surfnorm(x,y,z);

if disegna==1,
  surfl(x,y,z,'light');
  hold on
  shading interp;
  colormap('copper');
  daspect([1 1 1])
  camlight 
  lighting gouraud;
  axis equal
  xlabel('x_m');
  ylabel('y_m');
  zlabel('z_m');
  chois=1;
  if n>=11,
      display('EGT warning: For normal computation you have to wait a little bit!!');
      display('Press 1 to visulaize surface normals or 0 to abort');
      chois=input('(1/0)>> ');
  end;
  if chois==1,
      hold on
      if ogg==2,
          nx=-nx;
          ny=-ny;
          nz=-nz;
      end;
      quiver3(x,y,z,nx,ny,nz);
  end;
end;
