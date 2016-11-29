%%
% Converts a generic Conic to the Parametric form
%
%Input - inConic - input Conic
%           Conic form is [  A   B/2  D/2
%                           B/2   C   E/2
%                           D/2  E/2   F];
%
%Output - x,y - center of the ellipse
%         a - semiaxis
%         b - semiaxis
%         alpha_radians - tilt of the ellipse in radians
%
%%

function [x,y,a,b,alpha_radians] = f_conic2Param(inConic,toleranceForImConic)

size_inp = size(inConic);
if size_inp(1)==3 && size_inp(2)==3,
    % Normalization to ensure that the last element of the Conic is +/-1
%     Conic = inConic/abs(inConic(3,3));
    %Conic = inConic/inConic(1,1);
    A = inConic(1,1);
    B = inConic(1,2)*2;
    C = inConic(2,2);
    D = inConic(1,3)*2;
    E = inConic(2,3)*2;
    F = inConic(3,3);
else
    A = inConic(1);
    B = inConic(2);
    C = inConic(3);
    D = inConic(4);
    E = inConic(5);
    F = inConic(6);
end
%Converts the Conic vector of parameters to the Parametric form
[~,a,b,R,t]  = imconic([A B C D E F],toleranceForImConic);
 x = t(1);
 y = t(2);
 alpha_radians = atan2(R(2,1),R(2,2));
end

