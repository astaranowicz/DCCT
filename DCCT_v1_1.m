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
display('1: Load/Run Conf. Paper Data (limited-set)')
display('2: Load/Run User-input Data')
userInput = input('Select an option (default [any other key] exits): ');

run DCCT_variables_SetupScript
global DCCT_variables

if userInput == 1
    clc
    display('Running Conf. Paper Data (over 25 sphere images)')
    DCCT_calibration_method(1);
    
elseif userInput ==2
    clc
    display('User-input Data')
    display('  ')
    prompt = 'Enter an RGB sensor calibration matrix (e.g., [500 0 300; 0 500 250; 0 0 1]) (default - uses Dataset):  ';
    Kr = input(prompt);
    if ~isempty(Kr)
        DCCT_variables.Kr = Kr;
    end
    display('  ')
    display('DCCT User-input Options:')
    display('1: Load/Run Previously Saved Data (Loads Depth_userInput.mat and RGB_userInput.mat)')
    display('2: Input/Run New Data')
    display('  ')
    userInputOptions = input('Select an option(default [any other key] exits):  ');
    if userInputOptions == 1
        DCCT_calibration_method(2);
    elseif userInputOptions == 2
        prompt = 'Enter number of pictures to be used in brackets (e.g., [1 2 30]):  ';
        numOfImages = input(prompt);
        strprompt = 'Enter directory extension (default - /Data):  ';
        dirExt = input(strprompt, 's');
        if ~isempty(dirExt)
            DCCT_variables.dirExt = dirExt;
        end
        if ~isempty(numOfImages)
            DCCT_variables.setOfSpheres = numOfImages;
        end
        f_RGB_imageProcessing();
        f_Depth_imageProcessing();
        DCCT_calibration_method();
    end
else
    display('Exiting')
end

