%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v2.0 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% f_reflmatrix(Normal,Distance);
%
% Syntax:
% ------
%     Normal = plane normal vector
%     Distance = Distance from the mirror to the world reference system
%
%     H = reflection matrix 
%
% Description: 
% -----------
%     This function create the reflection matrix respect to plane dfined by
%     Normal and Distance
%       
% Example:
% -------   
%     close all; clear all
%     N = [1, 1, 0]; ds = 5;
%     H = f_reflmatrix(N,ds);  
%
% Author:
%    Stefano Scheggi
% Last update:
%    May, 2008
%
function H = f_reflmatrix(Normal,Distance);  
if nargin<2,
    display('EGT error: function "f_reflmatrix" needs 2 parameter at least');
elseif nargin>2,
    display('EGT warning: too much input parameters in "f_reflmatrix"!');
end;
I = eye(3);
H = double([(I-2*Normal*Normal') 2*Distance*Normal; 
               [0 0 0 ]   1  ]);


