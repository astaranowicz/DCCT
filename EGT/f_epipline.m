%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
%
% [la,ld] = f_epipline(Ua,Ud,F,flag,figa,figd,car);
%
% Syntax:
% ------
%      la(ld)= "epipolar lines (column vector)"
%      Ua(d) = "3xN matrix of N points in 3d scene of actual*(desired) camera"
%          F = "Fundamental matrix"
%       flag = "enables plot of epipolar lines"
%    figa(d) = "index of figure corresponding to the actual (desired) camera
%               frame"
%        car = "string containing color and shape of lines" 
%
% Descr: 
% ----- This function computes the pencil of epipolar lines. 
%       The pin-hole model is used. 
%
% Example:
% -------
%   clear all; close all;figure(2); hold on; figure(3); hold on;
% 	X=[0 , 5]; Y=[10, 4]; Z=[10,-3]; P=[X;Y;Z];     
% 	Rd=eye(3); td=[0,0,0]'; Hd=f_Rt2H(Rd,td);
% 	Ra=rotoy(-pi/6); ta=[-5,-5,0]'; Ha=f_Rt2H(Ra,ta);
% 	Kd=eye(3); Ka=eye(3);
% 	[ud,vd]=f_perspproj(P,Hd,Kd); 
% 	[ua,va]=f_perspproj(P,Ha,Ka);
% 	
% 	[ea,ed,F]=f_epipole(Ha,Hd,Ka,Kd);
% 	figure(2); grid on
% 	title('EGT- Epipolar Geometry - Actual Image plane and epipolar lines')
% 	plot(ea(1),ea(2),'rO'); text(ea(1)+.05,ea(2),'Epipole')
% 	plot(ua,va,'k*'); text(ua+.05,va,'Feature point')
% 	
% 	figure(3); grid on
% 	title('EGT- Epipolar Geometry - Desired Image plane and epipolar lines')
% 	plot(ed(1),ed(2),'gO'); text(ed(1)+.05,ed(2),'Epipole')
% 	plot(ud,vd,'k*'); text(ud+.05,vd,'Feature point')
%   Ua=[ua;va]; Ud=[ud;vd]; [la,ld]=f_epipline(Ua,Ud,F);
%
% Author:
%    Gian Luca Mariottini December 2005
function [la,ld]=f_epipline(Ua,Ud,F,flag,figa,figd,car);
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
     
   Oa([1:length(Ua(1,:))])=1; %estende con tanti "uni" quanta è la dim. di Ua
   ld=F*[Ua(1,:);Ua(2,:);Oa];
   
   Od([1:length(Ud(1,:))])=1; %estende con tanti "uni" quanta è la dim. di Ua
   la=F'*[Ud(1,:);Ud(2,:);Od];
   
   ea=null(F);
   ea=ea/ea(3);
   ed=null(F');
   ed=ed/ed(3);
   if flag==1,
       figure(figa);
       for i=1:length(la(1,:)),
            ma(i)=-la(1,i)/la(2,i);
            qa(i)=-la(3,i)/la(2,i);
            ya([1],i)=ma(i)*Ua(1,i)+qa(i);
            ya([2],i)=ma(i)*ea(1)+qa(i);
            plot([Ua(1,i) ea(1)],[ya(1,i) ya(2,i)],car);        
       end;
       figure(figd);
       for i=1:length(ld(1,:)),
            md(i)=-ld(1,i)/ld(2,i);
            qd(i)=-ld(3,i)/ld(2,i);
            yd([1],i)=md(i)*Ud(1,i)+qd(i);
            yd([2],i)=md(i)*ed(1)+qd(i);       
            plot([Ud(1,i) ed(1)],[yd(1,i) yd(2,i)],car);
       end;
   end;
   
   