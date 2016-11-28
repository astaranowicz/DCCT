%%
% Uses the Quadric cost function to calculate the Conic from center of the
% sphere
%
% Input - U_depth - set of points that belong on the sphere
%         Kd - depth calibration matrix
%         Dd - depth distortion vector
%         Kr - RGB calibration matrix
%         Dr - RGB distortion vector
%         R - Rotation matrix R_R_D
%         t - translation vector R_t_D
%
% Output - Conic - conic from the quadric cost function
%%


function Conic = f_Quadric2Conic(U_depth,Kd,Dd,Kr,Dr,R,t)

for i = 1:length(U_depth)
    %Constructing the 3D point vector to calculate the center
    X_depth(i).points = f_depth2XYZ(Kd,Dd,U_depth(i).points);
    M = f_sphereLinLS(X_depth(i).points(1:3,:));
    centerSphere_hat(i).center = M(1:3);
    radius_hat(i) = M(4);
end

%TODO: Utilize Dr
P = Kr*[R t];
for i = 1:length(U_depth)
    % Quadric
    Q = [eye(3) -centerSphere_hat(i).center;
        -centerSphere_hat(i).center'  centerSphere_hat(i).center'*centerSphere_hat(i).center-radius_hat(i)^2];
    
    Conic(i).conic = inv(P*inv(Q)*P');
end


end
