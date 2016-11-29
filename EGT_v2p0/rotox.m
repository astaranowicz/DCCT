function [m] = rotox(x)

    m=[1,0,0;0,cos(x),-sin(x);0,sin(x),cos(x)];