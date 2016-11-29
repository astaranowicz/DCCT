%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example by Gian Luca Mariottini
%
    close all
    clear all
    figure(1)
    title('Example 1 - 3D visualization of a thorus')
    grid on
    axis equal
    hold on
    f_3Dwf('k',3,'_{wf}');
    hold on
    view(-33,34);
    title('Epipolar Geometry Toolbox - 3D visualization - Example 1')
    hold on
    f_3Dsurface(2,rotoy(pi/6),[10,10,0]',[3 1 1],30,1); %torus