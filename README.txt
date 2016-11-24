README FILE


DCCT v.1.1 is a contribution that is part of the following paper:
Journal Version:
Aaron N. Staranowicz, Garrett R. Brown, Fabio Morbidi, Gian-Luca Mariottini, Practical and accurate calibration of RGB-D cameras using spheres, Computer Vision and Image Understanding, Available online 9 April 2015, ISSN 1077-3142, http://dx.doi.org/10.1016/j.cviu.2015.03.013.
Conference Version:
A. Staranowicz, G.R. Brown, F. Morbidi, and G.L. Mariottini. Easyto-Use and Accurate Calibration of RGB-D Cameras from Spheres. In R. Klette, M. Rivera, and S. Satoh, editors, Proc. 6th Pacific-Rim Symposium on Image and Video Technology, volume 8333, pages 265–278. Springer, 2014. 

%%%%%%%%%%%
Please note that the feature selection is not user-friendly.  We are updating and testing the latest verison which is described in the above paper.  

For testing purposes, this toolbox was used to generate all estimated calibration parameters used in the above paper for our method, however, there might be some bugs in this toolbox which could lead to inaccurate calibration parameters.  To ensure an accurate calibration, select the ellipses very carefully.

%%%%%%%%%%%

Authors’ Contact Information:
Aaron Staranowicz (Corresponding Author) and Gian-Luca Mariottini
email:{aaron.staranowicz@mavs.uta.edu, gianluca@uta.edu}
Dept. of Comp. Science and Eng., University of Texas at Arlington, Arlington, Texas


Please note, this zipped file (DCCT_v1.zip) contains parts of the code used to generate the figures and results of our calibration algorithm and code that allows users to input their own data to calibrate a depth camera.  A newer code version will include a GUI for better usability. 

This zipped file contains several toolboxes:
- A conic toolbox (imconic):
  --Copyright © 2010, Levente Hunyadi All rights reserved. NOTE:  Some changes were made to this toolbox but some of the underlying ideas of Levente Hunydai are still used.

- The Epipolar Geometry Toolbox  (EGT)  toolbox (v2.0):
  -- By: G.L. Mariottini and D. Prattichizzo

-A Hough Transform for Circles 
 -- Copyright ©, Yuan-Liang Tang, Dept. of Info. Management, Chaoyang University of Technology


The 4 main files are:
 DCCT_v1_1.m
f_RGB_imageProcessing.m
f_Depth_imageProcessing.m
DCCT_variables_SetupScript.m

The main script is called:  DCCT_v1_1.m   (which is the script to run in MATLAB)

The main script allows users to either run the data to obtain the calibration results and some of the figures used in the paper, or to input users’ own images of spheres along with corresponding depth maps.  Users will follow the options in Matlab’s command window.  

The user can either use previously saved data from a past session of using DCCT or start with a new set of data.  If the user wants to collect the data from a new set of images, 2 more options will be shown which are:  the number of images to be used and where they are located on the disk.  Then the two functions (f_RGB_imageProcessing.m and f_Depth_imageProcessing.m) will be called to collect the data needed for the calibration.

The function f_RGB_imageProcessing.m allows is used to process the RGB images obtained from a depth camera. 

The function works as follows:   The user will be asked to select a rectangle around the sphere (as close to the sphere’s perceived edge) and clicks inside the rectangle twice.  The cropped image will be processed by canny edge and hough transform for circles.  Usually only 1 circle is found, but in the case of multiple circles, the user will be asked to choose the best circle that will be processed.  The circle will have a band created around the edge pixels that belong to the sphere.  Those points are then used in an RANSAC-based ellipse fitting algorithm.  This estimated ellipse is then shown to the user and the user will either accept the estimated ellipse or start over with the same image.  The options will allow users to decrease or increase the threshold and sigma used by Canny Edge detection.  Note that in most cases the threshold needs to be decreased but the sigma can remain the default value.  The process starts again with the user selecting the rectangle around the sphere. 

NOTE:  A good estimated ellipse will as close as possible to the edge of the sphere.  A bad estimated ellipse would have a side far away from the edge of the sphere.

After all the images have been processed, the user is asked if they would like to save the data or not save the data.


The function f_Depth_imageProcessing.m allows is used to process the depth map obtained from a depth camera. 

The function works as follows:  The user will be asked to select a rectangle around the sphere’s perceived edge and then press enter on the keyboard.  The pixel points inside the rectangle will be processed to remove the pixel points that either belong to the background or have a depth of zero.  The user will then either accept the pixel points shown with an overlay of blue on the sphere, or reject and start again.  If the user rejects then a new maximum depth needs to be selected to either allow more or less pixel points to be selected.

NOTE:  A good selection of the pixel points would have the pixel points that belong to the sphere colored blue.  A bad selection would have selected pixel points that do not belong to the sphere and would belong to the floor or background.

DCCT_variables_SetupScript.m is a set of default values used throughout the calibration process.  These should not be changed.  If the user decides to change the default values then expect some or all parts of the code will not work and errors will be created.  


The data should be formatted in the following way:
RGB images should be named RGBImage_ball_#.jpg
Depth Maps should be named ball_#.txt (which includes the distance for each pixel in a column-wise vector which can be reshaped to match the RGB image.)
