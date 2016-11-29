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
% DCCT_toolbox v1.1
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

