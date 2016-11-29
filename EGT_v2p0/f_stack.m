function l=f_stack(M)

larg=length(M(1,:));
lung=length(M(:,1));
l=reshape(M,larg*lung,1);