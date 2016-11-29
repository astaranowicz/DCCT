%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Epipolar Geometry Toolbox v1.3 (EGT)    %
%                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  f_3Dpanoramic(H,color,quadric,a,b,r_rim)  
%
%USAGE:
%  Plot a 3D panoramic (parabolic/hyperbolic) camera; 
%  IMP= Note that the mirror frame is oriented as the Matlab frame
%
%SYNTAX:
%  H = Homogeneous matrix containing rotation R, and translation t, 
%      defined with respect to the world frame, of current camera. 
%      (for more details see the EGT manual).
%  color = character string which defines the colour used to plot the camera
%  quadric=1 ==> hyperbolic mirror
%         =2 ==> parabolic
%  a,b = hyperboloidal/parabolic mirror parameters in EGT frame according
%        to the following
%        HYPERBOLOID: (z+e)^2/a^2-(x^2+y^2)/b^2=1 where e=sqrt(a^2+b^2) 
%        PARABOLOID:  z+a/2= (x^2+y^2)/(2a)
%  r_rim = mirror parameter: it is the radius of the mirror rim
%
%EXAMPLE:
%    close all; clear all; figure(1); hold on; axis equal; view(58,18); 
%    Ha=[rotoy(0)*rotoz(0)*rotox(0) , [-0.2,-0.1,0]';
%               0    0    0            ,    1    ];
% 	 quadric=1; a=0.04; b=0.02; r_rim=0.04;  % Hyperbola
% 	 f_3Dwf('k',0.1); % World reference frame
% 	 f_3Dpanoramic(Ha,'g',quadric,a,b,0.03);   %Panoramic camera    
% 	 f_3Dframe(Ha,'g',0.06,'_{m}');   % Mirror reference frame
% 	 title('Example 1 - Panoramic Camera Placement and Plot'); grid on
%
%AUTHORS:
%  Gian Luca Mariottini- December 2005
%LAST UPDATE:
%  December 2005
%

function f_3Dpanoramic(H,color,quadric,a,b,r_rim);

if nargin==1,
    color = 'k';    
    quadric=1;
    a = 3;
    b = 1;
    r_rim = 2;
elseif nargin==2,
    quadric=1;
    a = 3;
    b = 1;
    r_rim = 2;
elseif nargin==3,
    a = 3;
    b = 1;
    r_rim = 2;
    fig = 1;
elseif nargin==4,
    b = 1;
    r_rim = 2;
elseif nargin==5,
    r_rim = 2;
elseif nargin>6,
    display('EGT error: too many inputs in f_3Dpanoramic')
end;
 
% H([1:3],[1:3]) = Rwf2egt*Regt2mir = Rwf2mir
Rrpy = H([1:3],[1:3]);   % Rrpy = Rwf2mir
t = H([1:3],4);  % t = twf2mir
       
%%%%%%  MIRROR DESIGN  %%%%%%
x = -(r_rim):.001:(r_rim); %Step for design the mirror: CHANGE IT if you need!
y = -(r_rim):.001:(r_rim);
[X,Y] = meshgrid(x,y);

if quadric==1,
    Z = (a/b)*sqrt(X.^2+Y.^2+b^2)-sqrt(a^2+b^2); % HYPERBOLIC MIRROR
    z_rim = (a/b)*sqrt(r_rim^2+b^2)-sqrt(a^2+b^2);
elseif quadric==2,
    Z=(X.^2+Y.^2)/(2*a)-a/2; % PARABOLIC MIRROR
    z_rim = (r_rim^2)/(2*a)-a/2;
end;



% from matrix X,Y,Z to vector Xvet,Yvet,Zvet
num = length(X(1,:));
for i=1:num,
     Xvet(num*(i-1)+1:i*num) = X(i,:);
     Yvet(num*(i-1)+1:i*num) = Y(i,:);
     Zvet(num*(i-1)+1:i*num) = Z(i,:);
end;   
   
EGTtrans = (Rrpy(:,1)*Xvet)+(Rrpy(:,2)*Yvet)+(Rrpy(:,3)*Zvet);
   for i=1:num,
      Xmatrix(i,:) = EGTtrans(1,num*(i-1)+1:num*i);  
      Ymatrix(i,:) = EGTtrans(2,num*(i-1)+1:num*i);  
      Zmatrix(i,:) = EGTtrans(3,num*(i-1)+1:num*i); 
   end   
   
% from EGT world frame to MATLAB frame
Rmat_egt = eye(3); 
MATtrans = (Rmat_egt(:,1)*EGTtrans(1,:))+(Rmat_egt(:,2)*EGTtrans(2,:))+(Rmat_egt(:,3)*EGTtrans(3,:));
   for i=1:num
      XmatrixMAT(i,:) = MATtrans(1,num*(i-1)+1:num*i)+t(1);  
      YmatrixMAT(i,:) = MATtrans(2,num*(i-1)+1:num*i)+t(2);
      ZmatrixMAT(i,:) = MATtrans(3,num*(i-1)+1:num*i)+t(3);
   end  

surfl(XmatrixMAT,YmatrixMAT,ZmatrixMAT);
shading('interp');
colormap('gray');

%Camera over the mirror
if quadric==1,
    e=sqrt(a^2+b^2);
    Rcam=eye(3);
    tcam=Rrpy*[0;0;-2*e];%This depends on hyperboloid formulas
elseif quadric==2,
    Rcam=eye(3);
    tcam=Rrpy*[0;0;-5*a]; %This is a fixed value and can be arbitrarily changed!
end    
scale=r_rim/3;
f_3Dcamera([  Rrpy          ,   t+tcam;  0  0   0  1  ],color,scale*2/3);
f_3Dframe([  Rrpy          ,   t+tcam; 0  0   0  1  ],color,4*scale);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MIRROR SUPPORT DESIGN
upvertex = [ -r_rim,  r_rim,  r_rim,  -r_rim ;
              r_rim,  r_rim, -r_rim,  -r_rim ;
              z_rim,  z_rim,  z_rim,   z_rim];

downvertex = [    -r_rim,     r_rim,       r_rim,       -r_rim    ;
                   r_rim,     r_rim,      -r_rim,       -r_rim    ;
              (z_rim+2*r_rim),  (z_rim+2*r_rim),   (z_rim+2*r_rim),      (z_rim+2*r_rim)];              


% from MIRROR frame to MATLAB frame
for i=1:length(upvertex(1,:))
    upvertex_mat(:,i)   = Rrpy*upvertex(:,i)+t;
    downvertex_mat(:,i) = Rrpy*downvertex(:,i)+t;
end;

up = [upvertex_mat, upvertex_mat(:,1)];
down = [downvertex_mat,   downvertex_mat(:,1)];
plot3(up(1,:),up(2,:),up(3,:),color)
plot3(down(1,:),down(2,:),down(3,:),color)
plot3([upvertex_mat(1,1),downvertex_mat(1,1)],[upvertex_mat(2,1),downvertex_mat(2,1)],[upvertex_mat(3,1),downvertex_mat(3,1)],color);
plot3([upvertex_mat(1,2),downvertex_mat(1,2)],[upvertex_mat(2,2),downvertex_mat(2,2)],[upvertex_mat(3,2),downvertex_mat(3,2)],color);
plot3([upvertex_mat(1,3),downvertex_mat(1,3)],[upvertex_mat(2,3),downvertex_mat(2,3)],[upvertex_mat(3,3),downvertex_mat(3,3)],color);
plot3([upvertex_mat(1,4),downvertex_mat(1,4)],[upvertex_mat(2,4),downvertex_mat(2,4)],[upvertex_mat(3,4),downvertex_mat(3,4)],color);
theta = 0:.1:2*pi;
circle = [r_rim*cos(theta); r_rim*sin(theta); z_rim*ones(1,length(theta))];
circleMAT = [Rrpy , t; 0 0 0 , 1]*[circle; ones(1,length(theta))];
plot3(circleMAT(1,:),circleMAT(2,:),circleMAT(3,:),color);