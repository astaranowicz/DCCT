%%
%Changes from the parameteric form to the foci of the ellipse
%
%
%Input - z1,z2 - center of the ellipse
%        a,b - major/minor axis length
%        alpha - tilt of the ellipse
%
%Output - F1 - focus of the first 
%         F2 - 2nd focus
%         c - distance between a and b
%%

function [F1,F2,c] = f_param2Foci(z1,z2,a,b,alpha)



c = sqrt(abs(a^2 - b^2));

F1 = [z1 - (cos(alpha)*c);
      z2 - (sin(alpha)*c)];

F2 = [z1 + (cos(alpha)*c);
      z2 + (sin(alpha)*c)];

end

