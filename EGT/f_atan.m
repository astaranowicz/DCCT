function u=f_atan(y,x)

if (y>=0),
    u=atan2(y,x);
    return
elseif y<0,
    u=2*pi+atan2(y,x); %+ because atan will be negative
    return
end;