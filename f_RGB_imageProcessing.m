%%
%This function processes the RGB Images to estimate the conic of the
%selected sphere.
%%

function f_RGB_imageProcessing()
global DCCT_variables

try
	load UserInput_RGB
catch
	%Do nothing
end

counter = 1;
cannyThreshold_croppedImage = 0.4;
cannySigma_croppedImage = 0.15;

cannyThreshold_Image = 0.15;
cannySigma_Image = 0.25;

redo = false;
show_circles = false;
%img_scale = 0.5;
enlargment_ratio = 0.3;

while counter ~= length(DCCT_variables.setOfSpheres)+1
    
    i = DCCT_variables.setOfSpheres(counter);
    imageName = [DCCT_variables.dataPath,DCCT_variables.RGBImageName,num2str(i),DCCT_variables.fileExt];
    image = imread(imageName);
    figure(1)
    imshow(image)
    hold on

	if redo || ~exist('Ellipse_RGB') || i > size(Ellipse_RGB, 2)
    	title(char(['RGB Image of sphere ', num2str(i)],['Select a box around the sphere profile'],['then press enter to continue']));
		circles_count = 0;
		while circles_count == 0
			if ~redo
				AA = imsubtract(image(:,:,1), rgb2gray(image));
				AA = medfilt2(AA, [3 3]);
				AA = im2bw(AA, 0.15);
				AA = bwareaopen(AA,300);
				AA = bwlabel(AA, 8);
				stats = regionprops(AA, 'BoundingBox', 'Centroid');
				if length(stats) > 0
					rect = stats(1).BoundingBox;
					%Make it a little bigger
					rect(1:2) = rect(1:2) - enlargment_ratio * rect(3:4);
					rect(3:4) = (1 + 2 * enlargment_ratio) * rect(3:4);
				end
			end
			if redo || length(stats) == 0
				rectHandle = imrect;
				title(char(['RGB Image of sphere ', num2str(i)],['If you finished selecting the box, press enter']));
				input('If you finished selecting the box, press enter','s');
				rect = rectHandle.getPosition;
			end
			[AA,rect]= imcrop(image, rect);
			rectangle('Position',rect,'EdgeColor', 'r')
			%AA = imresize(AA, img_scale);
			AA_bin = adaptivethreshold(AA, 15, 0.02, 1);
			display('Please wait while detecting the RGB circle...');
			circles = houghcircles(AA_bin,floor(0.35*size(AA_bin,1)), ceil(0.65*size(AA_bin,1)));
			if size(circles, 1) == 0
				AA_bin = rgb2gray(AA);
				circles = houghcircles(AA_bin,floor(size(AA_bin,1)/2*.70), ceil(size(AA_bin,1)/2*1.30));
			end
			circles_count = size(circles, 1);
			if circles_count == 0
				display('No circles were found! Try again...');
				redo = true;
			end
		end
		BW = edge(AA_bin,'canny',cannyThreshold_croppedImage,cannySigma_croppedImage);

		if show_circles
			figure(2);
			imshow(BW);
			axis image;
			hold on;
			title('Canny Edge Detection for Cropped Image')

			for k = 1: size(circles, 1),
				plot(circles(k,1), circles(k,2), 'r+');
				DrawCircle(circles(k,1), circles(k,2), circles(k,3), 32, 'b-');
			end
		end

		BW2 = edge(rgb2gray(image),'canny',cannyThreshold_Image,cannySigma_Image);

		%circles(:,1) = circles(:,1) / img_scale;
		%circles(:,2) = circles(:,2) / img_scale;
		%circles(:,3) = circles(:,3) / img_scale;
		%circles(:,4) = circles(:,4) / img_scale;

		circles(:,1) = circles(:,1)+rect(1);
		circles(:,2) = circles(:,2)+rect(2);


		if show_circles
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
		end

		if size(circles(:,1),1) == 1 || (~show_circles && size(circles(:,1),1) >= 1),
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

		Ellipse_RGB(i).t = [z1;z2];
		Ellipse_RGB(i).a = a;
		Ellipse_RGB(i).b = b;
		Ellipse_RGB(i).alpha = alpha;
    end
        
    projected_sphere_center_RGB(i).points = f_projectionSphere(Ellipse_RGB(i).t(1), Ellipse_RGB(i).t(2), Ellipse_RGB(i).a, Ellipse_RGB(i).b, Ellipse_RGB(i).alpha, DCCT_variables.Kr, DCCT_variables.ProjectSphereCentertolerance);
    
    
    fighandle = figure(1);
    imshow(image);
    hold on
    title(char(['Estimated Ellipse on top of Image ',num2str(i)], ['Is the current ellipse ok: y,n [y]']));
    Conic = f_param2Conic_Ellipse(Ellipse_RGB(i).t(1),Ellipse_RGB(i).t(2),Ellipse_RGB(i).a,Ellipse_RGB(i).b,Ellipse_RGB(i).alpha);
    f_drawConic(Conic,fighandle,'r');
    
    userResStr = lower(input('Is the current ellipse ok: y,n [y]:  ','s'));
    
    if ~isempty(strfind(userResStr,'n'))
        clc
        display(['Redoing image number ',num2str(counter)]);
        cannyThreshold_ImageResp = input('New Canny Threshold: (default: 0.15)  ');
        if isempty(cannyThreshold_ImageResp)
            cannyThreshold_Image = 0.15;
        else
            cannyThreshold_Image = cannyThreshold_ImageResp;
        end
        cannySigma_ImageResp = input('New Canny Sigma: (default: 0.25) ');
        if isempty(cannySigma_ImageResp)
            cannySigma_Image = 0.25;
        else
            cannySigma_Image = cannySigma_ImageResp;
        end
        houghCircles_Resp = lower(input('Show hough circles candicated to select: y,n [n]: ', 's'));
        show_circles = ~isempty(strfind(houghCircles_Resp,'y'));
        redo = true;
    else
        cannyThreshold_Image = 0.15;
        cannySigma_Image = 0.25;
        show_circles = false;
        redo = false;
        
        DCCT_variables.Ellipse_RGB(i) = Ellipse_RGB(i);
        DCCT_variables.projected_sphere_center_RGB(i) = projected_sphere_center_RGB(i);
        counter = counter + 1;
    end
    
    close all
end

userSaveRes = lower(input('Do you want to save this set: y,n [y]:  ','s'));
if ~isempty(strfind(userSaveRes,'n'))
	display('Not saving');
else
	display('Saving');
	save('UserInput_RGB','Ellipse_RGB','projected_sphere_center_RGB');
end

end







