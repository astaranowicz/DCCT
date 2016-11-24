%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
%
% function U=f_perspproj(X,H,K,adjplot);
%
% Syntax:
% ------
%    X = "matrix of points in 3d scene" (see manual) (can be also in homogeneous notation)
%    H = "Homogeneous matrix containing rotation R and translation t w.r.t. the world frame."
%    K = "Internal Parameters of the camera"
%    adjplot = 0 - "do not plot nothing"
%              1 - "plot of points in the image plane"
%              2 - "plot of lines connecting the camera center and the 3D
%                   feature (use only if you are visualizing a 3D camera setup)"
%
% Descr: 
% ----- This function computes the perspective projection of a set of feature points 
%       X expressed in the world-frame.
% 
% Ex:
% --
%     clear all
%     close all    
%     P=[0;10;10]; Kd=eye(3);     
%     figure(1); hold on; grid on; axis equal;  view(-45,14)   
%     f_3Dwf('k',3,'_{wf}'); f_scenepnt(P,'r*');             
%     Rd=rotoy(0); td=[0,-10,0]'; Hd=f_Rt2H(Rd,td);     
%     f_3Dframe(Hd,'b:',3,'_{c}'); f_3Dcamera(Hd,'b',1);     
%     Ud=f_perspproj(P,Hd,Kd,2);
%     title('Epipolar Geometry Toolbox - Example 3 - Perspective Projection');
%     figure(2); plot(ud,vd,'gO'); grid on;
%     title('Epipolar Geometry Toolbox - Example 3 - Image plane');
%
% Author:
%    Gian Luca Mariottini 
% Last Update:
%    July 2008
%      - two for cycles are joined for faster running
%
function U=f_perspproj(U,H,K,adjplot);

if nargin==3,
    plottatutto=0;
elseif nargin==4,
    plottatutto=adjplot;
end;

      R=H([1:3],[1:3]); %Matrice rotazione riferita alla rotaz. del robot su se stesso (roll-pitch-yaw)
      t=H([1:3],4); %vettore traslazione riferito al sistema wf.
      
      %Control U as homogeneous
      if length(U(:,1))==3, %not homogeneous
          Uo=[U;ones(1,length(U(1,:)))];
      else
          Uo=U;
      end
            
      Rw2i=R';
      tw2i=-Rw2i*t;
      
%       for i=1:length(Uo(1,:)),
%           projpnt(:,i)=K*[Rw2i tw2i]*[Uo([1:3],i) ; 1]; %Points must be firstly transformed in the egt frame
%           projpntnorm(1,i)=projpnt(1,i)/projpnt(3,i);
%           projpntnorm(2,i)=projpnt(2,i)/projpnt(3,i);
%           projpntnorm(3,i)=projpnt(3,i)/projpnt(3,i);
%             if plottatutto==2,
%                 Pc(:,i)=[Rw2i tw2i]*[Uo([1:3],i) ; 1];
%                 plot3([t(1) U(1,i)],[t(2) U(2,i)],[t(3) U(3,i)],'r:');  
%             end;  
%       end;

          P=K*[Rw2i tw2i];
          n=length(Uo(1,:));
          projpnt([1:3],:)=P*[Uo([1:3],:) ; ones(1,n)]; %Points must be firstly transformed in the egt frame
          u=projpnt(1,:)./projpnt(3,:);
          v=projpnt(2,:)./projpnt(3,:);
            if plottatutto==2,
                Pc([1:3],:)=[Rw2i tw2i]*[Uo([1:3],:) ; ones(1,n)];
                for i=1:n
                    plot3([t(1) U(1,i)],[t(2) U(2,i)],[t(3) U(3,i)],'r:');  
                end
            end;
            U=[u;
               v];
            
            
%       u=projpntnorm(1,:);
%       v=projpntnorm(2,:);
      
      % Plot of lines joining camera and 3D features
        
          
          