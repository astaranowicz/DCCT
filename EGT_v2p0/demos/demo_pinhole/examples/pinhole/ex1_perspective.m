%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.1 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
%   
% -* Perspective projection example *- 
%       Placement of camera and plot of point on the image plane
%
     clear all
     close all
     X=0;     
     Y=10;     
     Z=10;     
     P=[X;Y;Z];     
     figure(1);     
     hold on; grid on; 
     axis equal;     
     f_3Dwf('k',3,'_{wf}');     
     f_scenepnt(P,'r*');             
     Rd=rotoy(0)*rotox(0)*rotoz(0);       
     td=[0,-10,0]';     
     Hd=f_Rt2H(Rd,td);     
     f_3Dframe(Hd,'b:',3,'_{c}');     
     f_3Dcamera(Hd,'b',1.5,2);     
     Kd=eye(3);     
     Ud=f_perspproj(P,Hd,Kd,2);  % put "2" to see the line joining 
                                      % the camera center and the feature   
     title('Epipolar Geometry Toolbox - Example 3 - Perspective Projection');
     for i=1:2:360,
         view(i,24)
         pause(.1)
     end    