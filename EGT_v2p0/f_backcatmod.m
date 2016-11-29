function xs_b=f_backcatmod(m1,m,ell,K);
if length(m1(:,1))==2,
    for i=1:length(m1(1,:)),
        m1(3,i)=1;
    end;
end
for i=1:length(m1(1,:)),
    xt=inv(K)*m1(:,i);
    xt(3)=-(ell+m); %% Must be forced because the thrd component of m1 is 1 and remains 1 also in xt
    lambda=(ell*(ell+m) + sqrt(-(ell^2-1)*(xt(1)^2+xt(2)^2)+(ell+m)^2))/(xt(1)^2+xt(2)^2+(ell+m)^2);
    cp=[0;0;ell];
    xs_b(:,i)=lambda*xt+cp;
end