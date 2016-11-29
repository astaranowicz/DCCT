%%
%Converts the pixel points to the 3D points
%
% Input - Kd - depth camera calibration matrix
%         ud - pixel points with form: [u,v,Z]
%
% Output - XYZ - 3D point in the depth camera
%%

function XYZ = f_depth2XYZ(Kd,ud)

Pd = [inv(Kd) zeros(3,1);
      zeros(1,3)  1];
u = ud(1,:);
v = ud(2,:);
Z = ud(3,:);
%As in eq (6) from the paper
%form: [uZ,vZ,Z];
ud_bar = [u.*Z;v.*Z;Z;ones(1,length(ud(1,:)))];

XYZ = Pd *ud_bar;

end

