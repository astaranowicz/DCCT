function bw = adaptivethreshold(img,ws,C,tm)
%ADAPTIVETHRESHOLD An adaptive thresholding algorithm that seperates the
%foreground from the background with nonuniform illumination.
%  bw = adaptivethreshold(img,ws,C) outputs a binary image bw with the local 
%   threshold mean-C or median-C to the image img.
%  ws is the local window size.
%  tm is 0 or 1, a switch between mean and median. tm = 0 mean(default); tm = 1 median.
%
%  Contributed by Guanglei Xiong (xgl99@mails.tsinghua.edu.cn)
%  at Tsinghua University, Beijing, China.
%
%  For more information, please see
%  http://homepages.inf.ed.ac.uk/rbf/HIPR2/adpthrsh.htm

if (nargin<3)
    error('You must provide the image img, the window size ws, and C.');
elseif (nargin == 3)
    tm = 0;
elseif (tm~=0 && tm~=1)
    error('tm must be 0 or 1.');
end

if size(img,3) == 3
	img = rgb2gray(img);
end

if tm == 0
    mimg = imfilter(img,fspecial('average',ws),'replicate');
else
    mimg = medfilt2(img,[ws ws]);
end
simg = mimg-img-C;
bw = im2bw(simg,0);
bw = imcomplement(bw);
