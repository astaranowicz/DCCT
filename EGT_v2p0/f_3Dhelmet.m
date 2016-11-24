%% >> f_3Dhelmet(Hh_w,Hc_h,color,b_plot)
% 
% Description: This function plots an helmet made of "n" single pinhole
% cameras, all painted with a "color". "b_plot" is a boolean (1=plot, 2=no
% plot).
%  
% H = [Rh_w th_w; 
%         0'   1] *single camera case*
% is a matrix containing the homogeneous transf. matrix from the <wf> to
% the head frame <h>.
%
% In the case of multi-cameras we will have
%
% H = [H1]
%     [H2]
%      ...
%     [Hn]; %is a tensor!!
% where each Hi has the form of H described above.

function f_3Dhelmet(Hh_m, Hc_h,size,color,b_plot)
% determine the numebr of cameras present
if b_plot==1,
    n=(length(Hc_h(1,1,:)));
    for i=1:n,
       f_3Dcamera_helm(Hh_m, Hc_h(:,:,i), 0.1,color);
    end
end;