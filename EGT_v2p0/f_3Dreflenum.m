%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% f_2Dreflenum(Struct_Point, color, distance);
%
% Syntax:
% ------
%     Struct_Point = struct('point',[],'number',[]);
%     color = color of points
%     distance = distance from text and point    
%
% Description: 
% -----------
%     This function plot the points defined in Struct_Point.point
%     and their label defined in Struct_Point.number 
%
% Author:
%    Stefano Scheggi
% Last update:
%    Mar., 2008
%

function f_3Dreflenum(Struct_Point, color, distance);
    if nargin<1
        display('EGT error: function "f_3Dreflenum" needs 1 parameter at least');
    elseif nargin<2
        color = 'r*';
        distance = 0;
    elseif nargin<3
        distance = 0;
    elseif nargin>3
        display('EGT warning: too much input parameters in "f_3Dreflenum"!');
    end;

    for index = 1:size(Struct_Point, 2)
        point = Struct_Point(index).point;
        label = strcat(num2str(Struct_Point(index).number));
        plot3(point(1), point(2), point(3), color);
        text(point(1)+distance, point(2)+distance, point(3)+distance, label);
    end
