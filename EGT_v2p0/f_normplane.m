%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox  (EGT)  %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  f_normplane     Normal versor of a plane.
%
%  [plane,flag] = f_normplane(X1,X2);
%  X1,X2 = vectors belonging to the plane; the normal versor to the plane 
%          is defined as the cross product between X1 and X2.
%
%  Descr: 
%  -----  Given two vectors, this function computes the normal versor to the 
%         plane to which the vectors belong.
%  
%  plane = f_normplane(X1,X2) returns in the 3by1 vector 'plane' the normal
%  versor to the plane where the vectors X1 and X2 lie. If X1 and X2 are 
%  linearly dependent the plane is not defined, in this case the function 
%  assigns to the variable 'plane' the not valid value [0,0,0]'.
%  [plane,flag] = f_normplane(X1,X2) returns flag=0 if the two vectors 
%  X1 and X2 are linearly dependent or flag=1 if they are independent.
%  
%
%  Eleonora Alunno - January 2004


function [plane,flag] = f_normplane(X1,X2)

% control to verify if the points are homogeneous    
if length(X1(:,1))==3, % non omogeneo
     X1o = [X1;ones(1,length(X1(1,:)))];
  else
      X1o = X1;
      X1 = X1o(1:3,:);
end
    
if length(X2(:,1))==3, % non omogeneo
     X2o = [X2;ones(1,length(X2(1,:)))];
  else
      X2o = X2;
      X2 = X2o(1:3,:);
end

% test to verify if the vectors are linearly dependent or no
A = [ zeros(1,3) 1; X1o'; X2o' ];
% if rank(A,10^-08)<3
%     flag = 0;
%     plane = [0 0 0]';
%     disp('  EGT Warning: the following vectors are linearly dependent')
%     disp('               the plane is not defined');
%     display(X1);
%     display(X2);
% else
    flag = 1;
    plane = f_skew(X1)*X2;
%     plane = plane/norm(plane);
% end;
    %%