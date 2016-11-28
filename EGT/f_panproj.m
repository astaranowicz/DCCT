%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox v1.3 (EGT)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [q,Xhmir]=f_panproj(X,H,K,a,b,quadric,color)
%
%  Catadioptric projection of the scene points (hyperbolic/parabolic) mirror.
%
%  X = matrix containing the scene points, it can be a 4-by-n matrix 
%      (if the points are in homogeneous coordinates) or 3-by-n matrix 
%      (no homogeneous coordinates) where n is the number of points.  
%  H = the homogeneous transformation with respect to the world reference 
%      frame containing both rotation R and translation t, of current 
%      camera. H is a 4by4 matrix (for more details, see the EGT manual).
%  K = calibration matrix with the internal parameters of the camera. 
%  quadric = 1 for hyperbolic mirror;
%            2 for parabolic mirror.
%  a,b = mirror parameters:
%        HYPERBOLOID: (z+e)^2/a^2-(x^2+y^2)/b^2=1 where e=sqrt(a^2+b^2);
%        PARABOLOID:  z+a/2= (x^2+y^2)/(2a).
%  color = character string which defines the colour used to draw.
%
%  Descr: 
%  -----  This function projects the scene points in X into image points
%         through a panoramic camera with hyperbolic mirror and returns the
%         pixel coordinates of the image points.
%
%  q = f_panproj(X) or more simply f_panproj(X), returns into the
%  2-by-n matrix 'q', the pixel coordinates of the image points, 
%  in which are projected the scene points by the camera (for default 
%  the panoramic camera is assumed to have the same orientation of the 
%  EGT frame and to be positioned in the origin of the matlab frame,
%  mirror parameters default values are a = 3 cm, b = 1 cm and the 
%  calibration matrix K is equal to identity).
%  [q,Xhmir] = f_panproj(X) returns, in addition to the pixel coordinates 
%  of image points, the projection of the scene points into the mirror 
%  surface expressed in the "mirror coordinate system". 
%  [q,Xhmir] = f_panproj(X,H,K,a,b,fig,color), in addition to what 
%  has been said above, permits to plot the mirror projections into 
%  the mirror surface. To do this, in 'fig' has to be the number of figure 
%  where is the plot of the camera for which the mirror projection are 
%  calculated. 
%
% Example:
%   close all; clear all; figure(1); hold on; f_3Dwf('k',0.1);
% 	H=[rotoy(0)*rotoz(0)*rotox(-180*pi/180) , [-0.2,-0.1,0]';
%             0    0    0       ,    1    ];
% 	quadric=2; a=0.03; b=1;    r_rim=0.05; %Parabola
% 	f_3Dpanoramic(H,'g',quadric,a,b,r_rim);     
% 	f_3Dframe(H,'g',0.06,'_{m}');  % Mirror reference frame
% 	%Scene Point
% 	X=[ 0  -.2  .1; 0   .1  .2; .15   .3  .1];
% 	% Camera calibration matrix
% 	K=[10^3   0   320;  0   10^3  240; 0  0  1];
% 	% Point projection
% 	figure(1); axis equal; view(69,18);
% 	[q,Xh]=f_panproj(X,H,K,a,b,quadric,'b:');
% 	title('Example 2 - Central Catadioptric Imaging of 3D scene points');
% 	plot3( X(1,:), X(2,:), X(3,:),'r*'); grid on
% 	f_3Dwfenum(X,'k',0.01);
%
%Authors:
%     Eleonora Alunno 
%     Gian Luca Mariottini
%Last update:
%     Nov, 2005

function [q,Xhmir] = f_panproj(X,H,K,a,b,quadric,col)
% ho tolto r_rim, non serve a meno che non debba considerare il cono
% d'ombra al di sotto dello specchio

   if nargin==1
       H = [eye(3),[0 0 0]';
            0 0 0 , 1];
       K = eye(3);
       a = 3;
       b = 1;
       quadric=1;
       col = 'r.';
   elseif nargin==2
       K = eye(3);
       a = 3;
       b = 1;
       quadric=1;
       col = 'r.';
   elseif  nargin==3
       a = 3;
       b = 1;
       quadric=1;
       col = 'r.';
   elseif nargin==4
       b = 1;
       quadric=1;
       col = 'r.';
   elseif nargin==5
       quadric=1;
       col = 'r.';
   elseif nargin==6
       col = 'r.';
   elseif nargin>7
       display('  EGT error: too many inputs in f_panproj')
   end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  control to verify if X is homogeneous
  if length(X(:,1))==3, % non omogeneo
      Xo = [X ; ones(1,length(X(1,:)))];
  else
      Xo = X;
  end

% rotation and traslation between sdr camera and sdr mirror  
    Rrpy = H([1:3],[1:3]);  % Rrpy = Regt_mir
    t = H([1:3],4);         % t = tmat_mir

    tmir_mat=H([1:3],4);
    Rmir_mat=H([1:3],[1:3]);
    Rmat_mir=Rmir_mat';
    tmat_mir=-Rmir_mat'*(tmir_mat);
    Rmir_cam=eye(3);
    tmir_cam = [0 0 2*sqrt(a^2+b^2)]';

Xmir=[ Rmir_mat'  , -Rmir_mat'*tmir_mat;
        0  0  0   ,           1        ]*Xo;   
    
for i=1:length(Xmir(1,:)),
    Xhmir([1:3],i) = f_projectmir(a,b,quadric,Xmir(1:3,i));
    Xhmir(4,i) = 1; %to make it homogeneous
end;

if quadric==1, % HYPERBOLOID
    panprojpnt=K*[Rmir_cam tmir_cam]*Xhmir;
    for i=1:length(panprojpnt(1,:)),
        panprojpntnorm(1,i)=panprojpnt(1,i)/panprojpnt(3,i);
        panprojpntnorm(2,i)=panprojpnt(2,i)/panprojpnt(3,i);
        panprojpntnorm(3,i)=panprojpnt(3,i)/panprojpnt(3,i);
    end;
        q(1,:)=panprojpntnorm(1,:);
        q(2,:)=panprojpntnorm(2,:);
        q(3,:)=ones(1,length(panprojpnt(1,:)));
elseif quadric==2,%PARABOLOID
    Xhmir_3byn=Xhmir([1:3],:);
    q=K*[ [1 0 0; 0 1 0]*Xhmir_3byn;
          ones(1,length(Xhmir_3byn(1,:))) ];
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  from mirror frame to matlab frame   
Xhmat = [Rmir_mat    tmir_mat; 
            [0 0 0]     1    ]*Xhmir;  % i.e., Xhmat = Rmat_mir*Xhmir+tmat_mir

% PLOT
    for i=1:length(X(1,:)),
        plot3(X(1,i),X(2,i),X(3,i),strcat(col,'*'));
        plot3(Xhmat(1,i),Xhmat(2,i),Xhmat(3,i),strcat(col,'.'));
        if length(col)==2,
            col2=col(1);    
        elseif length(col)==1,
            col2=col;
        elseif length(col)>2,
            display('EGT error: color parameter too long. Can be only of length 2')
        end
        plot3([X(1,i) t(1)],[X(2,i) t(2)],[X(3,i) t(3)],strcat(col2,':'));
    end

