%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.1 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
%
% function v=f_perspline(U,H,K,adjplot,fig);
%
% function ell=f_perspline(PL,H,K,adjplot,fig);
%
% Syntax:
% ------
%   U       = "matrix of 2 points in 3d scene belonging to the line"
%   H       = "Homogeneous matrix containing rotation R and translation t w.r.t. the world frame (with f_Rt2H()"
%   K       = "Internal Parameters of the camera"
%   PL      = "6x1 vector representing the plucker coordinates of the line
%              (in the camera frame)"
%   adjplot = "boolean variable enabling plot of lines"
%   ell     = "line parameters"
% Descr: 
% ----- This function computes the perspective projection of a line (in 3D
%       space).
% 
%
% Author:
%    Gian Luca Mariottini 
% Last Update:
%    December 2004

function usc=f_perspline(U,H,K,adjplot,fig);

if nargin==3,
    adjplot=0,
    fig=1;
elseif nargin==4,
    fig=1,
end;

if length(U(:,1))==6,
    plucker=1;
elseif length(U(:,1))==3,
    plucker=0;
else    
    display('EGT error in calling f_perspline');
    return
end;    

%% 3) Normal of the interpretation plane 
    t=H(1:3,4);
    R=H([1:3],[1:3]);
    Rw2c=R';          
    tw2c=-Rw2c*t;
    scale=1/4;%1/K(2,2)*200;
    Rc2w=R;
    tc2w=t;
    
    P1= U([1:3],1);    
    P2= U([1:3],2);
    
    X1c=scale*(Rw2c*P1+tw2c); %3D point in the CAMERA frame      
    Xh1m=(Rc2w*X1c+tc2w); %3D point in MAtlab frame (De-rotation+aligning with EGT & transl.)
    
    X2c=scale*(Rw2c*P2+tw2c);
    Xh2m=(Rc2w*X2c+tc2w);
    
if adjplot==1, 
    plot3(Xh1m(1),Xh1m(2),Xh1m(3),'b.');
    plot3(Xh2m(1),Xh2m(2),Xh2m(3),'b.');
    plot3([Xh1m(1),Xh2m(1)],[Xh1m(2),Xh2m(2)],[Xh1m(3),Xh2m(3)],'b');
end;   
    %2.1) Vector joining 2 points and plot
            u=[(P2(1)-P1(1));(P2(2)-P1(2));(P2(3)-P1(3))];
            if adjplot==1,
                quiver3(P1(1),P1(2),P1(3),1/2*u(1),1/2*u(2),1/2*u(3),'r');
            end;    
    %2.2) Normal vector        
            n=f_skew(X1c([1:3]))*X2c([1:3])/norm(f_skew(X2c([1:3]))*X1c([1:3]));
            nm=Rc2w*n; %write n in the matlab frame        
            tm_w=H([1:3],4);
            if adjplot==1,
                quiver3(tm_w(1),tm_w(2),tm_w(3),2*nm(1),2*nm(2),2*nm(3),'b'); %normal plot
                text(tm_w(1)+2*nm(1),tm_w(2)+2*nm(2),tm_w(3)+2*nm(3),'n_{\pi}');
            end;    
    %2.3) Line equation (parameters)
            [u,v]=f_perspproj(U,H,K);
            Mell=[u(1) v(1) 1;
                  u(2) v(2) 1]; %Matrix of the line
            ell=null(Mell); %line parameters

        usc=[ell,Xh1m,Xh2m];