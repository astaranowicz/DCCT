%
%Uses the 6 point algorithm to estimate the K matrix and
%Rotation/Translation.
%
%Input - 3D world points that are non-planar, U Points from Image
%Output - Estimated Rotation, translation, and Calibration

function [R_hat,t_hat,K_hat] = f_sixPointAlgo(centerSphere_hat,U2)%Image_points)


World_points = [centerSphere_hat(1).center,centerSphere_hat(2).center,centerSphere_hat(3).center,centerSphere_hat(4).center,centerSphere_hat(5).center,centerSphere_hat(6).center];
World_points = [World_points;ones(1,size(World_points, 2))];
numOfPoints = length(World_points);

%Image_points = [U2(1).points,U2(2).points,U2(3).points,U2(4).points,U2(5).points,U2(6).points];

Image_points = [U2.points];

%calculates A matrix from Ax = 0
A = [];
for i = 1:numOfPoints,
    
    ai = [0, 0, 0, 0, (-1* World_points(:,i)'), (Image_points(2,i)*World_points(:,i)'); 
         (1*World_points(:,i)'), 0, 0, 0, 0,(Image_points(1,i)*-1*World_points(:,i)');
         (Image_points(2,i)*-1*World_points(:,i)'),  (Image_points(1,i)*World_points(:,i)'), 0, 0, 0, 0];

     A = [A; ai];
end

testing_rank = rank(A);

if testing_rank < 11
    R_hat = zeros(3,3);
    t_hat = zeros(3,1);
    K_hat = zeros(3,3);
    display('ill-conditioned');
    keyboard
end

[~,~,V]=svd(A);

x_hat = V(:,end);
%stacks the vector x in to the matrix P 
% W_P_C = [x_hat(1),x_hat(2),x_hat(3),x_hat(4);
%          x_hat(5),x_hat(6),x_hat(7),x_hat(8);
%          x_hat(9),x_hat(10),x_hat(11),x_hat(12)];

W_P_C = [ x_hat(1:4)' ; 
          x_hat(5:8)' ;
          x_hat(9:12)'];%OK
     
% calculates the translation
W_T_C = (-1*inv(W_P_C(1:3,1:3))) * W_P_C(1:3,4);

[Q,R] = qr(inv(W_P_C(1:3,1:3)));
%calculates the Rotation matrix
W_R_C = Q;
% calculates the intrinsic parameters
K = inv(R);
%normalizes to the last element of the K matrix
K = K/K(3,3);
% Scan columns of the estimated matrix "K" to adjust sign of "W_R_C"
% for i=1:3,
%     signK=find(K(:,i)<0);
%     if ~isempty(signK) & (K(signK,i)>0),
%         K(:,i)=-K(:,i);
%         W_R_C(i,:)=-W_R_C(i,:);
%     end
%     clear signK
% end
% % Final check on det(C_R_W_hat)...invert its sign if necessary
% if det(W_R_C) < 0,%1,
%     W_R_C=-W_R_C;
% end



changes=0;
for i=1:3,
    if K(i,i)<0;   
        K(i,:) = - K(i,:); % if 1/fu or 1/fv are negative, then change the sign of that row
        W_R_C(:,i) = - W_R_C(:,i); % and change the corresponding column in R
        changes = changes + 1;
    end
end
% Final change of sign, in case only one column of R was changed of sign.
if (mod(changes,2)~=0)||(det(W_R_C) <0), 
    W_R_C = -W_R_C;
end



K_hat = K;
R_hat = W_R_C'; 
t_hat = -R_hat*W_T_C; 
end


