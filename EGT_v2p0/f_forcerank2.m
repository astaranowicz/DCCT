function M=f_forcerank2(Mi);
if (size(Mi,1)~=3) | (size(Mi,2)~=3),
   disp('EGT Error: incorrect number of parameters')
else
   [U,S,V] = svd(Mi);
   S(3,3)=0;
   M = U*S*V';
end