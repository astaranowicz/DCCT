%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v2.0 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% pointR = f_reflpoint(Normal,Distance,Point,plot_flag,color);
%
% Syntax:
% ------
%     Normal = plane normal vector
%     Distance = Distance from the mirror to the world reference system
%     Point = points [X1 X2 ... Xn;
%                     Y1 Y2 ... Yn;
%                     Z1 Z2 ... Zn]
%     plot_flag = if '1' plot the points      
%
%     pointR = Point reflected 
%
% Description: 
% -----------
%     This function reflect points respect to a plane defined by Normal and
%     Distance and plot the reflected points
%       
% Example:
% -------   
%     close all; clear all
%     N = [1, 1, 0]; ds = 5;
%     X = [-10;30;0];
%     Xr = f_reflpoint(N,ds,X,1,'r*'); 
%
% Author:
%    Stefano Scheggi
%    Gian Luca Mariottini
% Last update:
%    May., 2008
%
function pointR = f_reflpoint(Normal,Distance,Point,plot,color); 
if nargin<3,
    display('EGT error: function "f_reflpoint" needs 3 parameter at least');
elseif nargin<4
    plot = 0;
    color = 'r*';
elseif nargin<5
    color = 'r*';
elseif nargin>5,
    display('EGT warning: too much input parameters in "f_reflpoint"!');
end;

D = f_reflmatrix(Normal,Distance);  

for index = 1:length(Point(1,:))
    pointR(:,index) = D*[Point([1:3],index);1];
end

% Draw points %
if plot == 1
   f_scenepnt(pointR,color,1); 
end

               
               
               
               
               
               