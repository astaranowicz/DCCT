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
%This function processes the RGB Images to estimate the conic of the
%selected sphere.
%%

function f_RGB_imageProcessing()
global DCCT_variables

counter = 1;
cannyThreshold_croppedImage = 0.4;
cannySigma_croppedImage = 0.25;

cannyThreshold_Image = 0.25;
cannySigma_Image = 0.25;


while counter ~= length(DCCT_variables.setOfSpheres)+1
    
    i = DCCT_variables.setOfSpheres(counter);
    imageName = [DCCT_variables.dirExt,DCCT_variables.RGBImageName,num2str(i),DCCT_variables.fileExt];
    image = imread(imageName);
    figure(1)
    imshow(image)
    hold on
    title(char(['RGB Image of sphere ', num2str(i)],[' Select a box around the sphere profile'],[' then click with right mouse button to select "crop"']));
    
    [AA,rect]= imcrop(image);
    rectangle('Position',rect,'EdgeColor', 'r')
    circles = houghcircles(AA,floor(size(AA,1)/2*.70), ceil(size(AA,1)/2*1.30));
    BW = edge(rgb2gray(AA),'canny',cannyThreshold_croppedImage,cannySigma_croppedImage);
    
    figure(2);
    imshow(BW);
    axis image;
    hold on;
    title('Canny Edge Detection for Cropped Image')
    
    for k = 1: size(circles, 1),
        plot(circles(k,1), circles(k,2), 'r+');
        DrawCircle(circles(k,1), circles(k,2), circles(k,3), 32, 'b-');
    end
    
    BW2 = edge(rgb2gray(image),'canny',cannyThreshold_Image,cannySigma_Image);
    
    circles(:,1) = circles(:,1)+rect(1);
    circles(:,2) = circles(:,2)+rect(2);
    
    for k = 1 : size(circles, 1),
        figure(k+100);
        imshow(BW2);
        axis image;
        hold on;
        title(['Canny Edge Detection for Image ',num2str(i),',  Circle Detected Number ', num2str(k)]);
        plot(circles(k,1), circles(k,2), 'r+');
        DrawCircle(circles(k,1), circles(k,2),circles(k,3), 32, 'b-');
        H = text(circles(k,1), circles(k,2), int2str(k));
        set(H,'Color',[1,1,1])
    end
    
    if size(circles(:,1),1) == 1,
        x = circles(1,1);
        y = circles(1,2);
        r = circles(1,3);
    elseif size(circles(:,1),1) > 1,
        userInput = input('Input Number: ');
        x = circles(userInput,1);
        y = circles(userInput,2);
        r = circles(userInput,3);
        clear userInput;
    else
        return
    end
    theta = 0 : (2 * pi / 32) : (2 * pi);
    pline_x = r * cos(theta) + x;
    pline_y = r * sin(theta) + y;
    
    pline_x_out = (r*1.10) * cos(theta) + x;
    pline_y_out = (r*1.10) * sin(theta) + y;
    
    pline_x_in = (r*.80) * cos(theta) + x;
    pline_y_in = (r*.80) * sin(theta) + y;
    
    mask_out= poly2mask(pline_x_out,pline_y_out,size(BW2,1),size(BW2,2));
    mask_in = poly2mask(pline_x_in,pline_y_in,size(BW2,1),size(BW2,2));
    
    img_mask = mask_out.*(1-mask_in);
    
    filter_BW = BW2.*img_mask;
    find_vec = find(filter_BW == 1);
    
    [indices_c, indices_r]= ind2sub([size(filter_BW,1),size(filter_BW,2)],find_vec);
    
    [z1,z2,a,b,alpha,phi_RGB,in_notused,o_notused] = f_ellipseFit_points2Ellipse([indices_r indices_c]', DCCT_variables.threshold_RGB);
    
    Ellipse_RGB(counter).t = [z1;z2];
    Ellipse_RGB(counter).a = a;
    Ellipse_RGB(counter).b = b;
    Ellipse_RGB(counter).alpha = alpha;
    
    projected_sphere_center_RGB(counter).points = f_projectionSphere(z1, z2, a, b, alpha, DCCT_variables.Kr, DCCT_variables.ProjectSphereCentertolerance);
    
    
    han = figure(1);
    imshow(image);
    hold on
    title(['Estimated Ellipse on top of Image ',num2str(i)]);
    Conic = f_param2Conic_Ellipse(z1,z2,a,b,alpha);
    f_drawConic(Conic,han,'r');
    
    userResStr = input('Is the current ellipse ok: y,n [y]:  ','s');
    
    if ~isempty(strfind(userResStr,'n'))
        clc
        display(['Redoing image number ',num2str(counter)]);
        cannyThreshold_ImageResp = input('New Canny Threshold: (default: 0.25)  ');
        cannySigma_ImageResp = input('New Canny Sigma: (default: 0.25) ');
        if isempty(cannySigma_ImageResp)
            cannySigma_Image = 0.25;
        else
            cannySigma_Image = cannySigma_ImageResp;
        end
        if isempty(cannyThreshold_ImageResp)
            cannyThreshold_Image = 0.25;
        else
            cannyThreshold_Image = cannyThreshold_ImageResp;
        end
    else
        cannyThreshold_Image = 0.25;
        cannySigma_Image = 0.25;
        
        DCCT_variables.Ellipse_RGB(counter) = Ellipse_RGB(counter);
        DCCT_variables.projected_sphere_center_RGB(counter) = projected_sphere_center_RGB(counter);
        counter = counter+1;
    end
    
    close all
end

userSaveRes = input('Do you want to save this set: y,n [n]:  ','s');
if ~isempty(strfind(userSaveRes,'y'))
    display('Saving');
    save('RGB_userInput','Ellipse_RGB','projected_sphere_center_RGB');
else
    display('Not saving');
end

end







