%%
%  This function processes the Depth Map to find the pixel points that
%  belong to the sphere and removes the background.
%%


function f_Depth_imageProcessing()

global DCCT_variables
counter = 1;

while counter ~= length(DCCT_variables.setOfSpheres)+1
    
    i = DCCT_variables.setOfSpheres(counter);
    depthMapName = [DCCT_variables.dataPath,DCCT_variables.DepthImageName,num2str(i),DCCT_variables.fileExt];
    imageRead = imread(depthMapName);
    %Load image and get points
    figure(1);
    imshow(imageRead,[]);
    hold on
    title(char(['Depth Image of sphere ', num2str(i)],['Select a box around the sphere profile'], ['then press enter to continue']))
    %input for the edges of the selection box
    rectHandle = imrect;
    title(char(['Depth Image of sphere ', num2str(i)],['If you finished selecting the box, press enter']))
    input('If you finished selecting the box, press enter','s');
    rect = rectHandle.getPosition;
    Txmin = floor(rect(1));
    Tymin = floor(rect(2));
    Txmax = floor(rect(1)+rect(3));
    Tymax = floor(rect(2)+rect(4));
    %plots the rectangle that were clicked
    rectangle('Position',rect,'EdgeColor', 'r')
    %input to find the correct depth that the sphere is at
    depthFinder=[floor(rect(1)+(rect(3)/2)) floor(rect(2)+(rect(4)/2))];
    
    xGLM=imageRead([Tymin:Tymax],[Txmin:Txmax]);
    depth_click = imageRead(depthFinder(2),depthFinder(1));
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
    
    figure(1);
    title(char(['Depth Image of sphere ', num2str(i)],['Is the current selection of points on the sphere ok: y,n [y]']))
    userResStr = lower(input('Is the current selection of points on the sphere ok: y,n [y]:  ','s'));
    
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
        
        DCCT_variables.avg_depthPixelPoints(:,i) = avg_depthPixelPoints(:,i);
        DCCT_variables.U_depth_clicked(i) = U_depth_clicked(i);
        counter = counter+1;
    end
    close all
    
end

userSaveRes = lower(input('Do you want to save this set: y,n [y]:  ','s'));
if ~isempty(strfind(userSaveRes,'n'))
	display('Not saving');
else
	display('Saving');
	save('UserInput_Depth','avg_depthPixelPoints','U_depth_clicked');
end

end
