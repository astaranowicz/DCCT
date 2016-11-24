function u=f_atan(y,x)

%if (y>=0),
%    u=atan2(y,x)
%elseif y<0,
%    u=2*pi*ones(1,length(y)) + atan2(y,x); %+ because atan will be negative
%end;

u_1=atan2(y,x);
u_2=2.0*pi*ones(1,length(y)) + atan2(y,x); %+ because atan will be negative
ind=find(y<=0);
u=u_1;
if ~isempty(ind),
  u(ind)=u_2(ind);
end