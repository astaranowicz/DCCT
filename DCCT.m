%%
% DCCT_toolbox v1.1
% SUBMITTED: "Depth-Camera Calibration Toolbox (DCCT): accurate, robust, and
% practical calibration of depth cameras"
% Authors: Aaron Staranowicz, Fabio Morbidi, and Gian-Luca Mariottini
% (WAFR) 2012
%
% Please see the README File
%
%%
addpath('EGT_v2p0/')
addpath('Data/')
addpath('conic/')

clear all
close all
clc

display('DCCT Options:')
display('1: Load/Run Sample Data (over 50 sphere images)')
display('2: Load/Run User-input Data')
userInput = input('Select an option (default [any other key] exits): ');

global DCCT_variables
run DCCT_variables_SetupScript

if userInput == 1
    clc
    display('Running Sample Data (over 50 sphere images)')
    DCCT_calibration_method(1);
    
elseif userInput ==2
    clc
    display('User-input Data')
    display('  ')
    prompt = 'Enter an RGB sensor calibration matrix (e.g., [525 0 319.5; 0 525 239.5; 0 0 1]) (default: uses Dataset):  ';
    Kr = input(prompt);
    if ~isempty(Kr)
        DCCT_variables.Kr = Kr;
    end
    display('  ')
    display('DCCT User-input Options:')
    display('1: Load/Run Previously Saved Data (Loads UserInput_Depth.mat and UserInput_RGB.mat)')
    display('2: Input/Run New Data from hard disk')
    display('3: Input/Run New Data from Kinect sensor')
    display('  ')
    userInputOptions = input('Select an option(default [any other key] exits):  ');
    if userInputOptions == 1
        DCCT_calibration_method(2);
    elseif userInputOptions == 2
        prompt = 'Enter number of pictures to be used in brackets (e.g., [1 2 30]):  ';
        while true
			numOfImages = input(prompt);
			if isempty(numOfImages)
				numOfImages = DCCT_variables.setOfSpheres;
				break;
			end
			if length(numOfImages) >= 3
				break;
			end
			display('Enter at least three elements!')
		end
        DCCT_variables.setOfSpheres = numOfImages;
        strprompt = ['Enter directory to load images from (default: ', DCCT_variables.dataPath, '):  '];
        dataPath = input(strprompt, 's');
        if ~isempty(dataPath)
            DCCT_variables.dataPath = dataPath;
        end
        f_RGB_imageProcessing();
        f_Depth_imageProcessing();
        DCCT_calibration_method();
    elseif userInputOptions == 3
        grabFrames();
        f_RGB_imageProcessing();
        f_Depth_imageProcessing();
        DCCT_calibration_method();
    end
else
    display('Exiting')
end

