%%
%Compute the Euler Angles
%Input - R - rotation matrix
%
%Output -
%          row - z axis
%          pitch - y axis
%          yaw - x axis
%%


function [roll,pitch,yaw] = f_rotationMat2EulerAngles(R)

%From Siciliano Book
%The known parameters for the Rotation Matrix
%pitch range = (-pi/2, pi/2)

roll_est = atan2(R(2,1), R(1,1));
pitch_est = atan2(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));
yaw_est = atan2(R(3,2), R(3,3));
%pitch range = (pi/2, 3pi/2)
roll2_est = atan2(-R(2,1), -R(1,1));
pitch2_est = atan2(-R(3,1), -sqrt(R(3,2)^2 + R(3,3)^2));
yaw2_est = atan2(-R(3,2), -R(3,3));

if(pitch_est > -pi/2 && pitch_est < pi/2),
    roll = roll_est;
    pitch = pitch_est;
    yaw = yaw_est;
else
    roll = roll2_est;
    pitch = pitch2_est;
    yaw = yaw2_est;
end


end

