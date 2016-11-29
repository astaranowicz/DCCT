%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ROUNDN  Rounds input at n-th power of 10
%
%  y = ROUNDN(x) 
%  y = ROUNDN(x,n) rounds the input data x at a power
%  of tens.  n=-3 rounds the input data to
%  the 10^(-3) (thousand) position.
%

function [x] = f_roundn(x,n)

if nargin == 0
    error('EGT error: Incorrect number of arguments')
elseif nargin == 1
    n = -2;
end

res  = 10 ^ (fix(-n));

%  Set the significant digits for the input data

x = round(x * res) / res;
