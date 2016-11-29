function f_plot_conicwparams(a,b, x0, y0, alpha, col)

t=[0:1:360]*pi/180;
N=length(t);

x = x0*ones(1,N) + (a*cos(t)*cos(alpha)) - (b*sin(t)*sin(alpha));
y = y0*ones(1,N) + (a*cos(t)*sin(alpha)) + (b*sin(t)*cos(alpha));

plot(x,y,col)