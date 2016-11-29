%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox  (EGT) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function f_grad
%
% Computation of GRAD_L as in "The Foundamental matrix: theory,
% alghorithms and stability analysis", Luong and Faugeras eq. 18 e 19
% Normalization and  denormalization is taken into account.
%
% Author:
%     Gian Luca Mariottini 
%     Pucci Augusto
%     Ernesto di Iorio
%
% Last update:
%     May 2004
% Revised from:
%   Phil Torr - The Structure and Motion Toolkit for MATLAB.
%   http://wwwcms.brookes.ac.uk/ philiptorr
%
function F = f_grad(U,Up,F0)
    % features Hartley's normalization
    
    [nRows nCols]=size(U);

    %% Riporta F0 nel vettore a0
    f0= f_stack(F0);
    % ([Ud(:,1);1])'*F0*([Ua(:,1);1])

    %% trova a tale che minimizzi il funzionale di costo
    options1  = optimset('LevenbergMarquardt','on','TolX',1e-6,'TolFun',1e-11,'Display','off'); %6 e 11
    global U Up nCols;
    a=lsqnonlin(@cost,f0,[],[],options1);
    %% calcola la F da l vettore a -> Note that the rank(2) constraint is
    %% here enforced automatically
    F=unstkc(a,3,3);
return

% Distanza di Sampson (equazione 18)
function c=cost(fundmat)
    global U Up nCols;
    fundmat
    	Fe=unstkc(fundmat,3,3);
        Fe2=Fe;%f_forcerank2(Fe);
        f2=f_stack(Fe2);
    Fe2
    display('fine')
        if rank(Fe2)==1,
            display('rank=1')
            null(Fe2)
            pause
            alpha1=e(1);
            beta1=e(2);
            
            Fp = [f2(1),    f2(4),    alpha1*f2(1)+beta1*f2(4);
                  f2(2),    f2(5),    alpha1*f2(2)+beta1*f2(5);
                  f2(3),    f2(6),    alpha1*f2(3)+beta1*f2(6)];
        elseif (rank(Fe2)>=2)
            Fe2=f_forcerank2(Fe);
            e=null(Fe2);
            e(3);
            e=-e/e(3);
            alpha1=e(1);
            beta1=e(2);
            Fp = [f2(1),    f2(4),    alpha1*f2(1)+beta1*f2(4);
                  f2(2),    f2(5),    alpha1*f2(2)+beta1*f2(5);
                  f2(3),    f2(6),    alpha1*f2(3)+beta1*f2(6)];
        end
        
    c=0;
    for i=1:nCols,
        ld= (Fp*[U(:,i);1]);
        la= (Fp'*[Up(:,i);1]);
        normF = (la(1))^2 + (la(2))^2 + (ld(1))^2 + (ld(2))^2
        
        c=c+([Up(:,i);1]'*Fp*[U(:,i);1])^2/(normF);
    end