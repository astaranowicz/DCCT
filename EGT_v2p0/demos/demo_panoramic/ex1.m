%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1 - Central Catadioptric Camera Placement and Plot  

close all; clear all; figure(1); hold on; axis equal; view(58,18); 
Ha=[rotoy(0)*rotoz(0)*rotox(0) , [-0.2,-0.1,0]';
        0    0    0            ,    1    ];
quadric=1; a=0.04; b=0.02; r_rim=0.04;  % Hyperbola
%quadric=2; a=0.03; b=1;    r_rim=0.05; % Parabola
f_3Dwf('k',0.1); % World reference frame
f_3Dpanoramic(Ha,'g',quadric,a,b,0.03);   %Panoramic camera    
f_3Dframe(Ha,'g',0.06,'_{m}');   % Mirror reference frame
title('Example 1 - Panoramic Camera Placement and Plot'); grid on