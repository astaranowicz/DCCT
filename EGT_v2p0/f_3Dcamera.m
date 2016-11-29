%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% f_3Dcamera(H,color,scale,version);
%
% Syntax:
% ------
%         H = homogeneous trasformation matrix obtained by the use of f_Rt2H.
%     color = a string containing the color of lines of frame.
%     scale = axes scaling factor (magnitude)
%     version = {1,2} selects to visualize the std. camera (1) or the
%               pyramidal (transparent) on (2);
%
% Description: 
% -----------
%     This function plots a 3D pin-hole camera in 3D matlab frame.
%       
% Example:
% -------   
%     close all; clear all
%     figure(1); grid on; axis equal; hold on; f_3Dwf('k',2);
%     Rd=rotoz(0)*rotoy(pi/4)*rotox(0); td=[-10,-10,0]';  
%     Hd=f_Rt2H(Rd,td); f_3Dframe(Hd,'g',4,'_{cam}'); view(20,26);    
%     f_3Dcamera(Hd,'g',1); %3D pin-hole camera
%     f_3Dcamera(Hd,'r',1,2);
%     title('Epipolar Geometry Toolbox - 3D camera visualization ')
%
% Author:
%    Gian Luca Mariottini
% Last update:
%    Mar., 2007
%

function f_3Dcamera(H,color,scale,version,name_img);
if nargin==0
    display('EGT error: function "f_3Dcamera" needs 1 parameter at least');
elseif nargin==1,
    color='r';
    scale=1;
    version=1;
    name_img=[];
elseif nargin==2,
    scale=1;
    version=1;
    name_img=[];
elseif nargin==3,
    version=1;
    name_img=[];
elseif nargin>4,
    display('EGT warning: too much input parameters in "f_3Dcamera"!')
end;

if version<2,
    R=H([1:3],[1:3]);
    t=H([1:3],4);
    trasl=t;

    Rwf2m=[1,0,0;
           0,1,0;
           0,0,1];
    %CAmera points: to be expressed in the camera frame;
    CAMup=scale*[-1,-1,  1, 1, 1.5,-1.5,-1, 1;
                 1, 1,  1, 1, 1.5, 1.5, 1, 1
                2,-2, -2, 2,   3,   3, 2, 2];

    Ri2w=R;
    for i=1:length(CAMup(1,:)),    
      CAMupTRASF(:,i)=Ri2w*(CAMup(:,i))+trasl; %disegna la telec. attuale orientata nel verso giusto
    end;

    CAMdwn=scale*[-1,-1,  1, 1, 1.5,-1.5,-1, 1;
                -1,-1, -1,-1,-1.5,-1.5,-1,-1;
                 2,-2, -2, 2,   3,   3, 2, 2];
               
    for i=1:length(CAMdwn(1,:)),    
      CAMdwnTRASF(:,i)=Ri2w*(CAMdwn(:,i))+trasl; %disegna la telec. attuale orientata nel verso giusto
    end;

   %Riporto i punti dal sistema EGT al Matlab frame
   for i=1:length(CAMupTRASF(1,:)),    
      CAMupTRASFm(:,i)=CAMupTRASF(:,i);
      CAMdwnTRASFm(:,i)=CAMdwnTRASF(:,i);
   end;

    plot3(CAMupTRASFm(1,:),CAMupTRASFm(2,:),CAMupTRASFm(3,:),color);
    plot3(CAMdwnTRASFm(1,:),CAMdwnTRASFm(2,:),CAMdwnTRASFm(3,:),color);
    plot3([CAMupTRASFm(1,1),CAMdwnTRASFm(1,1)],[CAMupTRASFm(2,1),CAMdwnTRASFm(2,1)],[CAMupTRASFm(3,1),CAMdwnTRASFm(3,1)],color);
    plot3([CAMupTRASFm(1,2),CAMdwnTRASFm(1,2)],[CAMupTRASFm(2,2),CAMdwnTRASFm(2,2)],[CAMupTRASFm(3,2),CAMdwnTRASFm(3,2)],color);
    plot3([CAMupTRASFm(1,3),CAMdwnTRASFm(1,3)],[CAMupTRASFm(2,3),CAMdwnTRASFm(2,3)],[CAMupTRASFm(3,3),CAMdwnTRASFm(3,3)],color);
    plot3([CAMupTRASFm(1,4),CAMdwnTRASFm(1,4)],[CAMupTRASFm(2,4),CAMdwnTRASFm(2,4)],[CAMupTRASFm(3,4),CAMdwnTRASFm(3,4)],color);
    plot3([CAMupTRASFm(1,5),CAMdwnTRASFm(1,5)],[CAMupTRASFm(2,5),CAMdwnTRASFm(2,5)],[CAMupTRASFm(3,5),CAMdwnTRASFm(3,5)],color);
    plot3([CAMupTRASFm(1,6),CAMdwnTRASFm(1,6)],[CAMupTRASFm(2,6),CAMdwnTRASFm(2,6)],[CAMupTRASFm(3,6),CAMdwnTRASFm(3,6)],color);
    
else

    P_c=[0,  0, 0;
         1, -1, 2;
         1,  1, 2;
        -1,  1, 2;
        -1, -1, 2]*scale;
    n=length(P_c(:,1));
    
    switch color
        case 'r'
            color2=[.3 0 0];
        case 'g'
            color2=[0 .3 0];
        case 'b' 
            color2=[0 0 .3];
        case 'k'
            color2=[.5 .5 .5];
        case 'w'
            color2=[0 0 0];
    end
    

% These are the cameras 3D points belonging to  the head, expressed in the
% head frame...
% for the plot I have to express them back in the world frame
    Rh_w = H([1:3],[1:3]); %this is given as a rotation from <w> to <h> in reality
    Rw_m=eye(3);%rotox(-pi/2);
    Rh_m = Rw_m*Rh_w;
    th_m = H([1:3],4); %this instead in centered in the <M> frame
    Th_m = th_m*ones(1,n);
    P=(Rh_m*P_c' + Th_m)';

%COmponents of the 1st patch
X1=[P(2,1) ; P(3,1) ; P(4,1) ; P(5,1)];
Y1=[P(2,2) ; P(3,2) ; P(4,2) ; P(5,2)];
Z1=[P(2,3) ; P(3,3) ; P(4,3) ; P(5,3)];
%COmponents of the 2nd patch
X2=[P(1,1) ; P(2,1) ; P(3,1) ; P(1,1)];
Y2=[P(1,2) ; P(2,2) ; P(3,2) ; P(1,2)];
Z2=[P(1,3) ; P(2,3) ; P(3,3) ; P(1,3)];
%COmponents of the 2nd patch
X3=[P(1,1) ; P(3,1) ; P(4,1) ; P(1,1)];
Y3=[P(1,2) ; P(3,2) ; P(4,2) ; P(1,2)];
Z3=[P(1,3) ; P(3,3) ; P(4,3) ; P(1,3)];
%COmponents of the 2nd patch
X4=[P(1,1) ; P(4,1) ; P(5,1) ; P(1,1)];
Y4=[P(1,2) ; P(4,2) ; P(5,2) ; P(1,2)];
Z4=[P(1,3) ; P(4,3) ; P(5,3) ; P(1,3)];
%COmponents of the 2nd patch
X5=[P(1,1) ; P(5,1) ; P(2,1) ; P(1,1)];
Y5=[P(1,2) ; P(5,2) ; P(2,2) ; P(1,2)];
Z5=[P(1,3) ; P(5,3) ; P(2,3) ; P(1,3)];

%color=[.2 .2 .2];
opaqueness=.15;
fill3(X1,Y1,Z1,[.9 .9 .9],'FaceAlpha',opaqueness*3);
%warp(X1,Y1,Z1,);
fill3(X2,Y2,Z2,color2+[.5 .5 .5],'FaceAlpha',opaqueness);
fill3(X3,Y3,Z3,color2+[.3 .3 .3],'FaceAlpha',opaqueness);
fill3(X4,Y4,Z4,color2+[.1 .1 .1],'FaceAlpha',opaqueness);
fill3(X5,Y5,Z5,color2,'FaceAlpha',opaqueness);
end;