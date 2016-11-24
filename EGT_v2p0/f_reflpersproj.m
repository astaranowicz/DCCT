%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v2.0 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% [PointsProjected, PointsOnMirror, StructPointsProj, StructPointsOnMir] =
%       f_reflpersproj(Normal,Distance,Hworldcamera,Kcam,Point_3D,Precision,MarginX,MarginY,MarginZ);
%
% Syntax:
% ------
%     Normal = plane normal vector
%     Distance = Distance from the mirror to the world reference system
%     Hworldcamera = homogeneus matrix between world system and camera system
%     Kcam = cmera calibration matrix
%
%     PointsProjected = 2D points projected on image plane of the camera
%                       [X1 X2 .... Xn;
%                        Y1 Y2 .... Yn]
%     PointsOnMirror = 3D points projected on mirror surface
%                       [X1 X2 .... Xn;
%                        Y1 Y2 .... Yn;
%                        Z1 Z2 .... Zn]
%     StructPointsProj = struct('point',[],'number',[]) with labeled projected points 
%     StructPointsOnMir = struct('point',[],'number',[]) with labeled
%                         points on mirror
%     Precision = precisione used to make numerical confrontation
%     MarginX = X-margin of the mirror
%     MarginY = Y-margin of the mirror
%     MarginZ = Z-margin of the mirror
%
% Description: 
% -----------
%     This function project the reflected points on the camera image plane
%     and on the mirror surface
%
% Author:
%    Stefano Scheggi
%    Gian Luca Mariottini
% Last update:
%    May, 08
%

function [PointsProjected, PointsOnMirror, StructPointsProj, StructPointsOnMir] = f_reflpersproj(Normal,Distance,Hworldcamera,Kcam,Point_3D,precision,marginX,marginY,marginZ);
    if nargin<5
        display('EGT error: function "f_3Dreflpersproj" needs 5 parameter at least');
    elseif nargin==6,
        % Precision Value
        precision = 0.0001;
    elseif nargin>9
        display('EGT warning: too much input parameters in "f_3Dreflpersproj"!');
    end;
   
    % Vectors for results
    PointsProjected = []; 
    PointsOnMirror = [];
    
    % Structure for results
    StructPointsProj = struct('point',[],'number',[]);
    StructPointsOnMir = struct('point',[],'number',[]);
    
    % Virtual camera
    Hv_w = f_reflcamera(Normal,Distance,Hworldcamera);
    Ovirt = Hv_w*[0 0 0 1]';
    
    % Project 3D Points on MIRROR
    index_POM = 1;    
    for index = 1:length(Point_3D(1,:))
        
        % direction
        t = Ovirt([1:3],1)-Point_3D([1:3],index);
        
        % Intersection between line and mirror
        lambdaline = (Distance-(Normal(1)*Ovirt(1)+Normal(2)*Ovirt(2)+Normal(3)*Ovirt(3)))/(Normal(1)*t(1)+Normal(2)*t(2)+Normal(3)*t(3));
        
        % Line from virtual camera to 3D point on MIRROR
        tline = Ovirt(1:3,:)+lambdaline*t;
        
        % Check mirror margin
        if nargin<7
           % mirror projection 
          	StructPointsOnMir(index_POM).point = tline;
            StructPointsOnMir(index_POM).number = index;
            PointsOnMirror(:,index_POM) = tline;
           % Points projection
            reflPoint_3D = f_reflpoint(Normal,Distance,Point_3D,0);
            [u,v]=f_perspproj(reflPoint_3D([1:3],index),Hworldcamera,Kcam);
            StructPointsProj(index_POM).point = [u;v];
            StructPointsProj(index_POM).number = index;
            PointsProjected(:,index_POM) = [u;v];
           % index increment 
            index_POM = index_POM + 1;
        elseif nargin==7
           if (tline(1)>= marginX(1)-precision) && (tline(1)< marginX(2)+precision)
                StructPointsOnMir(index_POM).point = tline;
                StructPointsOnMir(index_POM).number = index;
                PointsOnMirror(:,index_POM) = tline;
              % Points projection
                reflPoint_3D = f_reflpoint(Normal,Distance,Point_3D,0);
                [u,v]=f_perspproj(reflPoint_3D([1:3],index),Hworldcamera,Kcam);
                StructPointsProj(index_POM).point = [u;v];
                StructPointsProj(index_POM).number = index;
                PointsProjected(:,index_POM) = [u;v];
              % index increment 
                index_POM = index_POM + 1;
           end;
        elseif nargin==8
           if (tline(1)>= marginX(1)-precision) && (tline(1)< marginX(2)+precision) && (tline(2)>= marginY(1)-precision) && (tline(2)< marginY(2)+precision)
                StructPointsOnMir(index_POM).point = tline;
                StructPointsOnMir(index_POM).number = index;
                PointsOnMirror(:,index_POM) = tline;
              % Points projection
                reflPoint_3D = f_reflpoint(Normal,Distance,Point_3D,0);
                [u,v]=f_perspproj(reflPoint_3D([1:3],index),Hworldcamera,Kcam);
                StructPointsProj(index_POM).point = [u;v];
                StructPointsProj(index_POM).number = index;
                PointsProjected(:,index_POM) = [u;v];
              % index increment 
                index_POM = index_POM + 1;
           end;
        elseif nargin==9
           if (tline(1)>= marginX(1)-precision) && (tline(1)< marginX(2)+precision) && (tline(2)>= marginY(1)-precision) && (tline(2)< marginY(2)+precision) && (tline(3)>= marginZ(1)-precision) && (tline(3)< marginZ(2)+precision)
                StructPointsOnMir(index_POM).point = tline;
                StructPointsOnMir(index_POM).number = index;
                PointsOnMirror(:,index_POM) = tline;
              % Points projection
                reflPoint_3D = f_reflpoint(Normal,Distance,Point_3D,0);
                [u,v]=f_perspproj(reflPoint_3D([1:3],index),Hworldcamera,Kcam);
                StructPointsProj(index_POM).point = [u;v];
                StructPointsProj(index_POM).number = index;
                PointsProjected(:,index_POM) = [u;v];
              % index increment 
                index_POM = index_POM + 1;
           end    
        end
    end

      
    