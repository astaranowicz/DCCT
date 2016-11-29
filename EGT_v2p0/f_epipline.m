%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v2.0 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
%
% [l,lp] = f_epipline(U,Up,F,flag,figa,figd,car);
%
% Syntax:
% ------
%      l (lp)= "epipolar lines (column vector)"
%      U (Up)= "3xN matrix of N points in 3d scene of actual*(desired) camera"
%          F = "Fundamental matrix s.t   ([Up;1])^T*F*[U;1]=0  "
%       flag = "enables plot of epipolar lines"
%    figa(d) = "index of figure corresponding to the actual (desired) camera
%               frame"
%        car = "string containing color and shape of lines" 
%
% Descr: 
% ----- This function computes the pencil of epipolar lines. 
%       The pin-hole model is used. 
%
%
% Author:
%    Gian Luca Mariottini, May 2008
function [l,lp]=f_epipline(U,Up,F,flag,figa,figd,car);
 if nargin==3,
     flag=1;
     figa=2;
     figd=3;
     car='r';
 elseif nargin==4,
     figa=2;
     figd=3;
     car='r';
 elseif nargin==5,
     figd=3;
     car='r';
 elseif nargin==6,
     car='r';
 elseif nargin>7
     display('EGT error: greater number of inputs in f_epipline !')
 end;
     
   O([1:length(U(1,:))])=1; %estende con tanti "uni" quanta è la dim. di Ua
   lp=F*[U(1,:);U(2,:);O];
   
   Op([1:length(Up(1,:))])=1; %estende con tanti "uni" quanta è la dim. di Ua
   l=F'*[Up(1,:);Up(2,:);Op];
   
   e=null(F);
   e=e/e(3);
   ep=null(F');
   ep=ep/ep(3);
   if flag==1,
       figure(figa);
       for i=1:length(l(1,:)),    
            plot([U(1,i) e(1)],[U(2,i) e(2)],car);        
       end;
       figure(figd);
       for i=1:length(lp(1,:)),
            plot([Up(1,i) ep(1)],[Up(2,i) ep(2)],car);
       end;
   end;
   
   