%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox  (EGT) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Syntax:
%-------
%   U = f_unstacked(M,dimension)
%
%Description:
%-----------
%   The vector "M" will be unstacked in a (dimension x dimension) matrix "U"
%
%Author
%------
% Gian Luca Mariottini
function unstacked = f_unstacked(M,dimension)

if mod(length(M),dimension)==0,
    howmany=length(M)/dimension;
  for i=1:howmany;
    if i==1,
        unstacked=M(1:howmany);
    else
        unstacked=[unstacked ,  M(1+howmany*(i-1):howmany+howmany*(i-1))];
    end;
  end  
else
   display('EGT error: the dimension of the vector to unstack does not match with the desired matrix dimension');
end    