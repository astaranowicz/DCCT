%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% f_3Dwf(color,scale,axislabel);
%
% Syntax:
% ------
%     color = a string containing the color of lines of frame (e.g., 'b')
%     scale = is a scale factor for world frame <wf>-axes
%     axislabel = string containing characters with which the axis will be labelled.
%
% Description:
% -----------
%     This function plots the world reference frame. 
%     Note that this frame is coincident with the 3D Matlab reference frame.    
%        
% Example 1:
% ---------   
%     figure(1); hold on; axis equal; grid on; %common figure settings
%     f_3Dwf('r',1,'_{myaxis}');
%
% Example 2:
% ---------
%     figure(1); hold on; axis equal; grid on; %common figure settings
%     f_3Dwf; % standard settings will be used by EGT!
%
% Author:
%   Gian Luca Mariottini, University of Siena, ITALY
% Last update:
%   Dec - 05
%
function f_3Dwf(color,scale,axislabel);

if nargin==0,
    color='k';
    scale=1;
    axislabel='_{wf}';  
elseif nargin==1,
    scale=1;
    axislabel='_{wf}';  
elseif nargin==2,
    axislabel='_{wf}';  
elseif nargin>3,
    display('EGT error: greater number of inputs in f_3Dwf!')
end;
%Sistema di Riferimento Telecamera 
    Owf=scale*[0,0,0]';
    Xwf=scale*[1,0,0]';
    Ywf=scale*[0,1,0]';
    Zwf=scale*[0,0,1]';

%Plot del sistema di riferimento camera (Gian)
 plot3([Owf(1),Xwf(1)],[Owf(2),Xwf(2)],[Owf(3),Xwf(3)],color)
 plot3([Owf(1),Ywf(1)],[Owf(2),Ywf(2)],[Owf(3),Ywf(3)],color)
 plot3([Owf(1),Zwf(1)],[Owf(2),Zwf(2)],[Owf(3),Zwf(3)],color)
 col=strcat(color,'>');
 plot3(Xwf(1),Xwf(2),Xwf(3),col);
 col=strcat(color,'V');
 plot3(Ywf(1),Ywf(2),Ywf(3),col);
 col=strcat(color,'^');
 plot3(Zwf(1),Zwf(2),Zwf(3),col);

%%% Inizio modifica Fabio plot frecce %%%
%hold on
%quiver3(Owf(1),Owf(2),Owf(3),Xwf(1),0,0,color);
%quiver3(Owf(1),Owf(2),Owf(3),0,Ywf(2),0,color);
%quiver3(Owf(1),Owf(2),Owf(3),0,0,Zwf(3),color);
%%% Fine modifica Fabio plot frecce %%%

text(Xwf(1),Xwf(2),Xwf(3),strcat('X',axislabel))
text(Ywf(1),Ywf(2),Ywf(3),strcat('Y',axislabel))
text(Zwf(1),Zwf(2),Zwf(3),strcat('Z',axislabel))