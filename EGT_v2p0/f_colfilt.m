%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function f_colfilt(I,U)
%           
% Descr: 
% -----  This function performs the color filtering process. A simple
%        mean-variance filtering on RGB levels is adopted.
%        This function needs a previous f_colparam call in order to 
%        select mean and variance values of the color to filter.
%        The wide of the bell around the mean value is adjusted with 
%        a value sigma_{r,g,b} for each component.
%
% Syntax:
% ------  
%       I = the image (RGB matrix);
%       U = set of parameters that can be retireved with f_colparam
%           and U=[mr,vr,mg,vg,mb,vb,sigma_r,sigma_g,sigma_b];
%       where mr,mg,mb = mean values of the RGB components of the color to filter; 
%             vr,vb,vg = variance values of the RGB components of the color to filter; 
%
% Author:
%     Gian Luca Mariottini 
% Thanks to:
%     Fabio Morbidi
%     Francesco Giallombardo
%
% Last update:
%    Dec. - 2005
%
function Ifilt=f_colfilt(I,U);
mr=U(1);
mg=U(2);
mb=U(3);
vr=U(4);
vg=U(5);
vb=U(6);
sigma_r=U(7);
sigma_g=U(8);
sigma_b=U(9);
Ibw=(((I(:,:,1)>=(mr-sigma_r*sqrt(vr)) & I(:,:,1)<=(mr+sigma_r*sqrt(vr)))) & ((I(:,:,2)>=(mg-sigma_g*sqrt(vg)) & I(:,:,2)<=(mg+sigma_g*sqrt(vg)))) & ((I(:,:,3)>=(mb-sigma_b*sqrt(vb)) & I(:,:,3)<=(mb+sigma_b*sqrt(vb))))) ;
Ifilt=~(Ibw);

