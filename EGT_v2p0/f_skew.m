%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox  (EGT) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function f_skew(vett)
%
% Descr: 
% -----  Computes the skew matrix of a vector "vett".
%
function A=f_skew(vett)
   A=[    0,   -vett(3),   vett(2);
       vett(3),    0,     -vett(1);
      -vett(2), vett(1),     0];