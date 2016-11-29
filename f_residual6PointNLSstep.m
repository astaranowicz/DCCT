%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C) 2016 Aaron Staranowicz and Gian Luca Mariottini
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%
%Residual for the 6 point algorithm cost function LS
%
%%

function [inliers,Residual_temp] = f_residual6PointNLSstep(R,t,pixelCenter,centersOfSphere,Kr, Inliers_6pnt,threshold)


temp_KrRt = Kr*[R t];


X_depth = centersOfSphere;

for i = 1:length(centersOfSphere)

    X_depth_temp = temp_KrRt * [X_depth(:,i);1];

    U_depth = X_depth_temp / X_depth_temp(3);

    Residual_temp(i) = norm([pixelCenter(:,i);1] - U_depth);
end


indices =  Residual_temp <= threshold;

inliers = Inliers_6pnt(find(indices));



end



