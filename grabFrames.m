%%
%  This function grabs new RGB-D frames from the Kinect sensor.
%%


function grabFrames()

global DCCT_variables

strprompt = ['Enter directory to save images to (default: ', DCCT_variables.dataPath, '):  '];
dataPath = input(strprompt, 's');
if ~isempty(dataPath)
	DCCT_variables.dataPath = dataPath;
end

display('Please wait while camera setup...');
rgb_dev = videoinput('kinect', 1);
depth_dev = videoinput('kinect', 2);

rgb_dev.FramesPerTrigger = 1;
depth_dev.FramesPerTrigger = 1;

rgb_dev.TriggerRepeat = Inf;
depth_dev.TriggerRepeat = Inf;

triggerconfig([rgb_dev depth_dev], 'manual');

start([rgb_dev depth_dev]);

i = 1;
fighandle = 1;
set(gcf, 'CurrentCharacter', '@');
while true
   	trigger([rgb_dev depth_dev]);
    % Get the acquired frames and metadata.
    [imgColor, ts_color, metaData_Color] = getdata(rgb_dev);
    [imgDepth, ts_depth, metaData_Depth] = getdata(depth_dev);

	if ~ishandle(fighandle)
		fighandle = figure('Name', 'Press [s] to store frame, [e] to stop capture', 'NumberTitle','off');
		set(gcf, 'Position', get(0,'Screensize'));
		subplot(1, 2, 1); title('RGB');
		subplot(1, 2, 2); title('Depth');
		hold on
	end
	subplot(1, 2, 1); image(imgColor);
	subplot(1, 2, 2); image(uint8(imgDepth / 32));

	shall_store = lower(get(gcf, 'CurrentCharacter'));
	if shall_store ~= '@'
		figure(fighandle);
		set(gcf, 'CurrentCharacter', '@');
	end
	if shall_store == 's'
		display(['Storying frame ', num2str(i), ' ...']);
		imwrite(imgColor, [DCCT_variables.dataPath, '/', DCCT_variables.RGBImageName, num2str(i), DCCT_variables.RGBFileExt]);
		imwrite(imgDepth, [DCCT_variables.dataPath, '/', DCCT_variables.DepthImageName, num2str(i), DCCT_variables.depthFileExt]);

		i = i + 1;
	elseif shall_store == 'e'
		if i > 3
			clc;
			commandwindow;
			userStopRes = lower(input('Are you sure you want to stop capture and start caliberation: y,n [n]:  ','s'));
			if ~isempty(strfind(userStopRes,'y'))
				break;
			end
		else
			warning('You need to capture at least three frames...');
		end
    end
end
close all;
delete([rgb_dev depth_dev]);
DCCT_variables.setOfSpheres = [1:i-1];

end
