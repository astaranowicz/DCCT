%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C) 2016 Aaron Staranowicz and Gian Luca Mariottini
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%  This function processes the Depth Map to find the pixel points that
%  belong to the sphere and removes the background.
%%


function f_Depth_imageProcessing()

global DCCT_variables
counter = 1;

while counter ~= length(DCCT_variables.setOfSpheres)+1
    
    i = DCCT_variables.setOfSpheres(counter);
    depthMapName = [DCCT_variables.dirExt,DCCT_variables.DepthMapName,num2str(i),DCCT_variables.fileExtDM];
    DepthMap = load(depthMapName);
    d = DepthMap(:,3); %in millimeters
    %Reshapes the vector of distances in to a matrix
    x3=reshape(d,640,480)';
    %Load image and get points
    figure
    imshow(x3,[]);
    hold on
    title(char(['Depth Image of sphere ', num2str(i)],['Select a box around the sphere profile'], ['then press enter to continue']))
    %input for the edges of the selection box
    rectHandle = imrect;
    pause
    rect = rectHandle.getPosition;
    Txmin = floor(rect(1));
    Tymin = floor(rect(2));
    Txmax = floor(rect(1)+rect(3));
    Tymax = floor(rect(2)+rect(4));
    %plots the rectangle that were clicked
    rectangle('Position',rect,'EdgeColor', 'r')
    %input to find the correct depth that the sphere is at
    depthFinder=[floor(rect(1)+(rect(3)/2)) floor(rect(2)+(rect(4)/2))];
    
    xGLM=x3([Tymin:Tymax],[Txmin:Txmax]);
    depth_click = x3(depthFinder(2),depthFinder(1));
    %Minimal distance that we can remove from the image
    Zmin = DCCT_variables.minDistanceFromCamera;
    %The center of the sphere plus some distance (in mm)
    Zmax = depth_click+DCCT_variables.maxDistanceFromCamera;
    [vGLM, uGLM] = find(xGLM>=Zmin & xGLM<=Zmax);
    zGLM = [];
    for j=1:length(vGLM),
        zGLM(j) = xGLM(vGLM(j),uGLM(j));
    end
    uGLM = uGLM + Txmin -1; %due to the image starts at 1 not 0
    vGLM = vGLM + Tymin -1; %due to the image starts at 1 not 0
    %plots the resulting points on the sphere
    plot(uGLM, vGLM, 'b.')
    %Stacks in to a 3 by N matrix, with last element converted from
    %millimeters to meters
    Points_G = [uGLM';vGLM';zGLM/1000];
    avgXdepth = mean(Points_G(1,:));
    avgYdepth = mean(Points_G(2,:));
    avgdepth = mean(Points_G(3,:));
    
    avg_depthPixelPoints(:,i) = [avgXdepth; avgYdepth; avgdepth];
    U_depth_clicked(i).points = Points_G;
    
    userResStr = input('Is the current selection of points on the sphere ok: y,n [y]:  ','s');
    
    if ~isempty(strfind(userResStr,'n'))
        clc
        display(['Redoing image number ',num2str(counter)]);
        depth_ImageResp = input('New Maximum Distance Threshold: (default: 300 millimeters)  ');
        
        if isempty(depth_ImageResp)
            maxDistanceFromCamera = 0.25;
        else
            maxDistanceFromCamera = depth_ImageResp;
        end
    else
        maxDistanceFromCamera = 0.25;
        
        DCCT_variables.avg_depthPixelPoints(:,counter) = avg_depthPixelPoints(:,counter);
        DCCT_variables.U_depth_clicked(counter) = U_depth_clicked(counter);
        counter = counter+1;
    end
    close all
    
end

userSaveRes = input('Do you want to save this set: y,n [n]:  ','s');
if ~isempty(strfind(userSaveRes,'y'))
    display('Saving');
    save('Depth_userInput','avg_depthPixelPoints','U_depth_clicked');
else
    display('Not saving');
end


end
