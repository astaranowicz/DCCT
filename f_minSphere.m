%%
%Cost function for a sphere for Non-Linear Least Squares
%
% Input - X0 - [x0,y0,z0,r] - estimated center of the sphere
%         PointsForNLS - [x,y,z] - the points around the sphere
%
%Output - F - vector to minimize with respect to x
%         J - Jacobian of F with respect to x
%
%%

function [F,J] = f_minSphere(X0)

 global PointsForNLS
 
 M = size(PointsForNLS,2);
 
 F = [];
 J = [];
%%
%Cost function
 G = PointsForNLS - X0(1:3)*ones(1,M);
 a = 9.8722e-6;%[m]
%  a0=125.0622;
%  a1=1.2666;
%  a2=334.4474;
 alpha_theta=1;
 alpha_phi=1;
 for i=1:M
     % Compute cartesian form of covariance matrix
     [R_car, Rpol, T] = f_RcarRpol(-PointsForNLS(:,i), a, alpha_theta, alpha_phi); %a0, a1, a2, alpha_theta, alpha_phi);
     sigma2_x =R_car(1,1);
     sigma2_y =R_car(2,2);
     sigma2_z =R_car(3,3);
     sigma2_w = ((PointsForNLS(:,i)-[X0(1); X0(2); X0(3)]).^2)'*[4*sigma2_x; 4*sigma2_y; 4*sigma2_z];
          
     F(i) = (norm(G(:,i))^2 - X0(4)^2)*(1/sqrt(sigma2_w));
 end
 
 
%%
%Building the Jacobian for a Sphere
 a = X0(1) - PointsForNLS(1,:);  %X component
 b = X0(2) - PointsForNLS(2,:);  %Y component
 c = X0(3) - PointsForNLS(3,:);  %Z component
 c_sq = sqrt(a.*a + b.*b + c.*c);  % The square of the 3 components

 
 J = [(a./c_sq)' (b./c_sq)' (c./c_sq)' -ones(M,1)];
 

end
