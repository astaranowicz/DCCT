%%
%
%%

function u_d = f_DepthImage2Depthmap(Image)


% From Image "I" to "u_d = [u;v;z]";     
n_r = length(Image(:,1));
n_c = length(Image(1,:));
U = combvec([1:n_c],[1:n_r]);   
It = Image';
Z = It(:);
u_d = [U;Z'];


end



