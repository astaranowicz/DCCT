%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Epipolar Geometry Toolbox  (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example by Gian Luca Mariottini

  clear all
  close all
  
  % Structured Scene (a cube)
  X = [-5, 5,  5, -5,-5, 5, 5 , -5];     
  Z = [ 5, 5, 15, 15, 5, 5,15 , 15];     
  Y = [15,15, 15, 15,25,25, 25, 25];     
  P = [X;Y;Z];      
  
  figure(1); title('Example 2 - Cube and random scene creation with EGT v1.3')
  hold on; axis equal;    
  f_3Dwf('b',3,'_{wf}');    
  f_scenepnt(P,'r*',1); 
  
  % Unstructered Scene : Random Points  
  P2=f_3Drandpoint(25,[10 10 10],4);
  f_scenepnt(P2,'r.',1); 
  
  grid on; view(-29,26)
 
  f_3Dwfenum(P,'k',.5);
  f_3Dwfenum(P2,'k',.5);
  
  