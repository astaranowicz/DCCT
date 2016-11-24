function R=f_Procrustes(X,Y)

%for i=1:length(X(1,:)),
%    M(:,:,i)=Y(:,i)*X(:,i)';%*Y(:,i)';
%end
%M_sum=sum(M,3);
M_sum=(Y')'*X';
[U,D,V]=svd(M_sum);
R=V*U';