%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.3 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% Sol = f_E2Motion(H)
%
% Syntax:
% ------
%     H = Essential matrix   
%     Sol = structure containing the 4 possible solutions of the decomposition   
%
% Description: 
% -----------
%     This function evaluates the 4 possible solution from the Essential matrix 
%
% Author:
%    Stefano Scheggi
% Last update:
%    Mar., 2008
%

function Sol = f_E2Motion(H)
    [u,s,v] = svd(H);
    sigma = (s(1,1)+s(2,2))/2;
    s = diag([sigma sigma 0]);
    
    %Controllare che u,v sono appartenenti a SO(3)
    
    t1_sk = u*rotoz(pi/2)*s*u';
    t1 = [t1_sk(3,2) -t1_sk(3,1) t1_sk(2,1)]';
    t2_sk = u*rotoz(-pi/2)*s*u';
    t2 = [t2_sk(3,2) -t2_sk(3,1) t2_sk(2,1)]';
    
    R1 = u*rotoz(pi/2)'*v';
    R2 = u*rotoz(-pi/2)'*v';
    
    if det(R1) < 0
        R1 = -R1;
    end
    if det(R2) < 0
        R2 = -R2;
    end
    
    Sol(:,:,1) = [R1, t1];
    Sol(:,:,2) = [R1, t2];
    Sol(:,:,3) = [R2, t1];
    Sol(:,:,4) = [R2, t2];