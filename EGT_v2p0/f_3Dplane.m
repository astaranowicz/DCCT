%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v2.0 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% eq_plane = f_3Dplane(Normal,Distance,xlimit,ylimit,zlimit);
%
% Syntax:
% ------
%     Normal = plane normal (column) vector 
%     Distance = Distance from the mirror to the world reference system (a
%                positive distance means that the plane is shifted of that distance
%                along the direction of the normal, otherwise is shifted in the
%                opposite direction of the normal)
%     xlimit = X-axis drawing range for PLOT
%     ylimit = Y-axis drawing range for PLOT
%     zlimit = Z-axis drawing range for PLOT
%
%     eq_plane = plane equation (symbolic) 
%
% Description: 
% -----------
%     This function plots a 3D plane in 3D matlab frame and 
%     gives the plane equation (sym) as output.
%       
% Example:
% -------   
%     close all; clear all
%     figure(1); grid on; axis equal; hold on; view(20,26);
%     N = [1, 1, 0]'; ds = 5;
%     eq_plane = f_3Dplane(N,ds,[-20,10],[-12,34],[2,5]);
%     title('Epipolar Geometry Toolbox - 3D plane visualization ')
%
% Author:
%    Stefano Scheggi
%    Gian Luca Mariottini
% Last update:
%    May, 2008
%
function eq_plane = f_3Dplane(Normal,Distance,xlimit,ylimit,zlimit);  
if nargin<2,
    display('EGT error: function "f_3Dplane" needs 2 parameter at least');
elseif nargin<3
    xlimit = [-10,10];
    ylimit = [-10,10];
    zlimit = [-10,10];
elseif nargin<4
    ylimit = [-10,10];
    zlimit = [-10,10];
elseif nargin<5
    zlimit = [-10,10];
elseif nargin>5,
    display('EGT warning: too much input parameters in "f_3Dplane"!');
end;

N=Normal/norm(Normal);
syms x y z;
P = [x y z].';
eq_plane = (N'*P)-Distance;
opaqueness=.15;

% plane function & drawing function
if(N(1)~=0)
    %xplane = solve(eq_plane, x);
    xplane = (Distance - N(2)*y - N(3)*z)/N(1);
    if(N(2)~=0)||(N(3)~=0)  
        [y,z] = meshgrid(ylimit, zlimit);
        surf(eval(xplane),y,z,'FaceAlpha',opaqueness); 
        return;
    else
        [y,z] = meshgrid(ylimit, zlimit);
        x=eval(xplane)*ones(length(z),length(y));
        surf(x,y,z,'FaceAlpha',opaqueness); 
        return;
    end;
elseif(N(2)~=0)
    yplane = (Distance - N(1)*x - N(3)*z)/N(2);
    if(N(3)~=0)
        [x,z] = meshgrid(xlimit, zlimit);
        surf(x,eval(yplane),z,'FaceAlpha',opaqueness); 
        return;
    else
        [x,z] = meshgrid(xlimit, zlimit);
        y=(Distance/N(2))*ones(length(x),length(z));
        surf(x,y,z,'FaceAlpha',opaqueness);  
        return;
    end;
elseif(N(3)~=0)
    zplane = (Distance - N(1)*x - N(2)*y)/N(3);
    [x,y] = meshgrid(xlimit, ylimit);
    z=(Distance/N(3))*ones(length(x),length(y));
    surf(x,y,z,'FaceAlpha',opaqueness); 
    return;
end; 