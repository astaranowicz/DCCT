%%
% Draws a conic using imconic's drawing function
%
% Input - Conic - the conic that will be drawn
%         fighandle - handle of the figure to draw the ellipse
%         color -  color of the ellipse edge
%%

function f_drawConic(Conic,fighandle, color)
A = Conic(1,1);
B = Conic(1,2)*2;
C = Conic(2,2);
D = Conic(1,3)*2;
E = Conic(2,3)*2;
F = Conic(3,3);

imconic([A B C D E F],fighandle,color,[]);

end


