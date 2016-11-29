  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox v1.1 (EGT) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% H = function f_Rt2H(R,t)
%
% R   is a 3x3 roll-pitch-yaw rotation matrix of camera frame wrt the EGT frame (see manual) 
% t   is a 3x1 translation vector
% H   is a 4x4 homogeneous transform. relating the camera frame with the world frame.
%
% Descr: 
% -----  This function generates the homogeneous matrix H corresponding to
%        the rotational matrix Rc2w and a translational vector tc2w.
%        The input traslational component t is the vector, centered in the
%        world frame, and pointing toward the camera frame. 
%        The input R is the rotation  matrix to bring the EGT frame on the 
%        camera frame (Rc2egt) (see manual) 
%
% Author:
% ------
%    Gian Luca Mariottini 
% Last update:
% -----------
%    December 2003
function H=f_Rt2H(R,t);
   if (size(R)~=[3 3])|(size(t)~=[3 1]),
       display('EGT error: input error in f_Rt2H. Please control dimensions.');
   end;
   H=[rotox(-pi/2)*R, t;    %NOTE: The rotational component 
      zeros(1,3)    , 1];