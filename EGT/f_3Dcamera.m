%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% f_3Dcamera(H,color,scale);
%
% Syntax:
% ------
%         H = homogeneous trasformation matrix obtained by the use of f_Rt2H.
%     color = a string containing the color of lines of frame.
%     scale = axes scaling factor (magnitude)
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
%     title('Epipolar Geometry Toolbox - 3D camera visualization ')
%
% Author:
%    Gian Luca Mariottini
% Last update:
%    Dec., 2005
%
function f_3Dcamera(H,color,scale);
if nargin==0
    display('EGT error: function "f_3Dcamera" needs 1 parameter at least');
elseif nargin==1,
    color='r';
    scale=1;
elseif nargin==2,
    scale=1;
elseif nargin>3,
    display('EGT warning: too much input parameters in "f_3Dcamera"!')
end;

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
plot3([CAMupTRASFm(1,2),CAMdwnTRASFm(1,2)],[CAMupTRASFm(2,2),CAMdwnTRASFm(2,2)],[CAMupTRASFm(3,2),CAMdwnTRASFm(3,2)],color)
plot3([CAMupTRASFm(1,3),CAMdwnTRASFm(1,3)],[CAMupTRASFm(2,3),CAMdwnTRASFm(2,3)],[CAMupTRASFm(3,3),CAMdwnTRASFm(3,3)],color)
plot3([CAMupTRASFm(1,4),CAMdwnTRASFm(1,4)],[CAMupTRASFm(2,4),CAMdwnTRASFm(2,4)],[CAMupTRASFm(3,4),CAMdwnTRASFm(3,4)],color)
plot3([CAMupTRASFm(1,5),CAMdwnTRASFm(1,5)],[CAMupTRASFm(2,5),CAMdwnTRASFm(2,5)],[CAMupTRASFm(3,5),CAMdwnTRASFm(3,5)],color)
plot3([CAMupTRASFm(1,6),CAMdwnTRASFm(1,6)],[CAMupTRASFm(2,6),CAMdwnTRASFm(2,6)],[CAMupTRASFm(3,6),CAMdwnTRASFm(3,6)],color)