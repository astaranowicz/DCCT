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



