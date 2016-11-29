%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v2.0 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% Hvirtual = f_reflcamera(Normal,Distance,Hworldcamera,plot,color,length,scale,title);
%
% Syntax:
% ------
%     Input:
%     Normal = plane normal vector
%     Distance = Distance from the mirror to the world reference system
%     Hworldcamera = homogeneus matrix between world system and camera system
%     plot = if '1' plot the points     
%     color = color of the camera
%     length = length of the camera frame
%     scale = scale factor of the camera
%       
%     Output:
%     Hvirtual = homogeneus matrix between world system and virtual camera system
%
% Description: 
% -----------
%     This function reflect the real camera respect to a plane defined by Normal and
%     Distance and plot the virtual camera
%       
% Example:
% -------   
%     close all; clear all
%     Twc = [15 18 0]';
%     Rwc = rotoy((-100*pi)/180);
%     Hwc=f_Rt2H(Rwc,Twc);
%     N = [1, 1, 0]; ds = 5;
%     Hwv = f_reflcamera(N,ds,Hwc,1,'g',6,10,'_{cam V}');
%
% Author:
%    Stefano Scheggi
% Last update:
%    May., 2008
%
function Hvirtual = f_reflcamera(Normal,Distance,Hworldcamera,plot,color,length,scale,title);
if nargin<3,
    display('EGT error: function "f_reflcamera" needs 3 parameter at least');
elseif nargin<4
    plot = 0;
    color = 'r';
    lenght = 0;
    scale = 1;
    title = '';
elseif nargin<6
    lenght = 0;
    scale = 1;
    title = '';
elseif nargin<7
    scale = 1;
    title = '';
elseif nargin<8
    title = '';
elseif nargin>8,
    display('EGT warning: too much input parameters in "f_reflcamera"!');
end;

D = f_reflmatrix(Normal,Distance);  
Hvirtual = D*Hworldcamera;

% Draw points %
if plot == 1
   f_3Dcamera(Hvirtual,color,scale);
   f_3Dframe(Hvirtual,color,length,title);
end