%%
%Cost function for minimizing the ellipse in non-linear least squares
%Geometric Fit in Parametric form:
%x = z + Q(alpha)x_prime
%x_prime = [a*cos(phi)
%           b*sin(phi)]
%
%Q(alpha) = [ cos(alpha) -sin(alpha)
%             sin(alpha) cos(alpha)]
%
%NonLinear Least Squares Problem
%g = [x1_i; x_2i] - [z1; z2] - Q(alpha)*[a*cos(phi_i); b*sin(phi_i)]
%
%
%
% Input - X0( described below)
%
% Output - F - the cost function to be minimized
%          J -  the jacobian of the cost function
%%

function [F,J] = f_minEllipse(X0)

global PointsForNLSEllipse

M = size(PointsForNLSEllipse,2);

phi = X0(1:M); %the angle for each point on the ellipse
alpha = X0(M+1); % the angle of the ellipse from the major axis
a = X0(M+2); % length of the ellipse on the x axis
b = X0(M+3); % length of the ellipse on the y axis
z1 = X0(M+4); % x component of center
z2 = X0(M+5); % y component of center

F = [];
J = [];

Q = [cos(alpha) -sin(alpha);
     sin(alpha)  cos(alpha)];
%Derivative of Q
Q_dot = [-sin(alpha) -cos(alpha);
          cos(alpha) -sin(alpha)];
 
%% From Gander: geometric fit in parametric form

r1_2 = PointsForNLSEllipse(1,:) - z1 - ((a*cos(phi)*cos(alpha))- (b*sin(phi)*sin(alpha)));
r2_2 = PointsForNLSEllipse(2,:) - z2 - ((a*cos(phi)*sin(alpha))+ (b*sin(phi)*cos(alpha)));   


F = [r1_2;
     r2_2];

%% Jacobian
   
G1_temp2 = -Q*[-a*sin(phi);b*cos(phi)];    
G2_temp2 = -Q_dot*[a*cos(phi);b*sin(phi)];
G3_temp2 = -Q*[cos(phi);zeros(1,M)];
G4_temp2 = -Q*[zeros(1,M);sin(phi)];

%Building the Original J
S2 = diag(G1_temp2(1,:));
C2 = diag(G1_temp2(2,:));

G12 = [S2 G2_temp2(1,:)' G3_temp2(1,:)' G4_temp2(1,:)' -ones(M,1) zeros(M,1)];
G22 = [C2 G2_temp2(2,:)' G3_temp2(2,:)' G4_temp2(2,:)' zeros(M,1) -ones(M,1)];

J_temp2= [G12;G22];
indexes2 = [1:2:2*M 2:2:2*M];
indexes22 = 1:2*M;
J(indexes2,:) = J_temp2(indexes22,:);


end
