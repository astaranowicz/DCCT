%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% f_3Dpointplanedistance(Normal,ds,Point);
%
% Syntax:
% ------
%     Normal = plane normal vector
%     ds = Distance from the mirror to the world reference system
%     Point = point in world coordinates
%
%     distance = distance from point to mirror
%
% Description: 
% -----------
%     This function find the distance between a point and a plane 
%       
% Example:
% -------   
%     close all; clear all
%     N = [1, 1, 0]; ds = 5;
%     P = [10,12,1];
%     distance = f_3Dpointplanedistance(N,ds,P);  
%
% Author:
%    Stefano Scheggi
% Last update:
%    Mar., 2008
%
function distance = f_3Dpointplanedistance(Normal,ds,Point);  
if nargin<3,
    display('EGT error: function "f_3Dpointplanedistance" needs 3 parameter at least');
elseif nargin>3,
    display('EGT warning: too much input parameters in "f_3Dpointplanedistance"!');
end;

N = Normal/norm(Normal);
distance = abs(N'*(N*ds-Point));