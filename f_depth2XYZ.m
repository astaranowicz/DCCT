%%
%Converts the pixel points to the 3D points
%
% Input - Kd - depth camera calibration matrix
%         Dd - depth camera distortion vector
%         ud - pixel points with form: [u,v,Z]
%
% Output - XYZ - 3D point in the depth camera
%%

function XYZ = f_depth2XYZ(Kd,Dd,ud)

Pd = [inv(Kd) zeros(3,1);
      zeros(1,3)  1];
u = ud(1,:);
v = ud(2,:);
Z = ud(3,:);
%As in eq (6) from the paper
%form: [uZ,vZ,Z];
ud_bar = [u.*Z;v.*Z;Z;ones(1,length(ud(1,:)))];

XYZ = Pd * ud_bar;

%Estimate distortion vector (by Hamdi Sahloul)
%TODO: Revise the correctness
r2 = XYZ(1,:).^2 + XYZ(2,:).^2;
dr = 1 + Dd(1) * r2 + Dd(2) * r2.^2 + Dd(5) * r2.^3;
dtx = 2 * Dd(3) * XYZ(1,:) .* XYZ(2,:) + Dd(4) * (r2 + 2 * XYZ(1,:).^2);
dty = 2 * Dd(4) * XYZ(1,:) .* XYZ(2,:) + Dd(3) * (r2 + 2 * XYZ(2,:).^2);
XYZ(1,:) = dr .* XYZ(1,:) + dtx;
XYZ(2,:) = dr .* XYZ(2,:) + dty;

end

