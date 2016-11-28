%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% f_3Dframe(H,color,scale,axislabel);
%
% Syntax:
% ------
%         H = homogeneous trasformation matrix obtained by the use of
%             f_Rt2H.
%     color = a string containing the color of lines of frame.
%     scale = define the size of axes (standard value=1)
% axislabel = string containing characters which the axis will be labelled.
%
% Descr: 
% -----  
%     This function plots the camera frame in 3D Matlab frame.    
%        
% Ex1:
% ---   
%     close all; clear all
%     figure(1); hold on; grid on;  axis equal
%     f_3Dwf('k');
%     Rd=rotoz(0)*rotoy(pi/4)*rotox(0); %Rd is the rotation wrt the EGT frame
%     td=[-4,-4,0]';  
%     Hd=f_Rt2H(Rd,td); %Hd is the homog. transf. wrt the World Frame
%     f_3Dframe(Hd,'g',1,'_{cam}'); %this is the camera frame
%     view(-33,34);
%     title('Epipolar Geometry Toolbox - 3D visualization - Example 1')
%     
%
% Author:
%   Gian Luca Mariottini
% Last update:
%   Nov, 29th -2005

function f_3Dframe(H,color,scale,axislabel);

R=H([1:3],[1:3]);% = Rw2i
t=H([1:3],4);

if nargin==1,
    color='r';
    scale=1;
    axislabel='_{c}';    
elseif nargin==2,
    scale=1;
    axislabel='_{c}';    
elseif nargin==3,
    axislabel='_{c}';    
elseif nargin>4,
    display('EGT error: greater number of inputs in f_3Dframe!')
end;

%Camera reference frame 
Oc=scale*[0,0,0]';
Xc=scale*[1,0,0]';
Yc=scale*[0,1,0]';
Zc=scale*[0,0,1]';

Ri2w=R;
Oc1=Ri2w*Oc+t;
Xc1=Ri2w*Xc+t;
Yc1=Ri2w*Yc+t;
Zc1=Ri2w*Zc+t;

xlabel('Xm');
ylabel('Ym');
zlabel('Zm');

plot3([Oc1(1),Xc1(1)],[Oc1(2),Xc1(2)],[Oc1(3),Xc1(3)],color)
plot3([Oc1(1),Yc1(1)],[Oc1(2),Yc1(2)],[Oc1(3),Yc1(3)],color)
plot3([Oc1(1),Zc1(1)],[Oc1(2),Zc1(2)],[Oc1(3),Zc1(3)],color)
col=strcat(color,'>');
plot3(Xc1(1),Xc1(2),Xc1(3),col);
col=strcat(color,'V');
plot3(Yc1(1),Yc1(2),Yc1(3),col);
col=strcat(color,'^');
plot3(Zc1(1),Zc1(2),Zc1(3),col);
text(Xc1(1),Xc1(2),Xc1(3),strcat('X',axislabel))
text(Yc1(1),Yc1(2),Yc1(3),strcat('Y',axislabel))
text(Zc1(1),Zc1(2),Zc1(3),strcat('Z',axislabel))

