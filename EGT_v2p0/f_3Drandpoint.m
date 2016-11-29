%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
%
% f_3Drandpoint(n,center,radius);
%
% Description:
% -----------
%    This function generates "n" spatial feature points lying inside a sphere of radius
%    known ("radius") and centered in a point [x y z] (i.e., the 1x3 vector "center").
% 
% Example:
% --------
%    P=f_3Drandnpoint(8,[0 20 5],7); % 8 feature points centered in [0 20
%                                    % 5] the radius of the spere is 7;
% Last update:
%   December 05
% Author:
%   Gian Luca Mariottini
%
function P=f_3Drandpoint(n,center,radius);
if nargin~=3,
    display('EGT error: the input arguments for "f_3Drandpoint" must be 3')
end;
for i=1:n,
    P(1:3,i)=[ center(1) + randn(1)*radius/3; center(2) + randn(1)*radius/3; center(3) + randn(1)*radius/3];
end;         