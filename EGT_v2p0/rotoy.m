function [m] = rotoy(x)

    m=[cos(x),0,sin(x);0,1,0;-sin(x),0,cos(x)];