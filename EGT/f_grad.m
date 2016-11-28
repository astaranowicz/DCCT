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
function F = f_grad(Ua,Ud)
    % features Hartley's normalization
    [Ua,Ta]=f_normalize(Ua);
	[Ud,Td]=f_normalize(Ud);
    [nRows nCols]=size(Ua);
    % Calcola il valore iniziale per la matrice Fondamentale
    F0=f_Festim(Ua,Ud,2);
    % Riporta F0 nel vettore a0
    % la scalatura può portare a divisioni per zero quindi è commentata
    % F0=F0/F0(3,3);
    a0(1)=F0(1,1);
    a0(2)=F0(1,2);
    a0(3)=F0(1,3);
    a0(4)=F0(2,1);
    a0(5)=F0(2,2);
    a0(6)=F0(2,3);
    astim=[a0(1) a0(4);a0(2) a0(5);a0(3) a0(6)]\(F0(3,:)');
    a0(7)=astim(1);
    a0(8)=astim(2);
    % trova a tale che minimizzi il funzionale di costo
    a=fminsearch(@cost,a0,optimset('LevenbergMarquardt','on'),Ua,Ud,nCols);
    % calcola la F da l vettore a -> Note that the rank(2) constraint is
    % here enforced automatically
    F=L2F(a);
    % denormalizzazione della matrice
	F=Td'*F*Ta;
return

% Distanza di Sampson (equazione 18)
function c=cost(a,Ua,Ud,nCols)
    Ft=L2F(a);
    c=0;
    for i=1:nCols
        c=c+([Ud(:,i);1]'*Ft*[Ua(:,i);1])^2/([1 1 0]*((Ft*[Ua(:,i);1]).^2) + [1 1 0]*((Ft'*[Ud(:,i);1]).^2));
    end
return

% Calcola la F dal vettore a (equazione 19)
function F=L2F(a)
    F=[ a(1)                        a(2)                 a(3) ;
        a(4)                        a(5)                 a(6) ;
        a(7)*a(1)+a(8)*a(4) a(7)*a(2)+a(8)*a(5) a(7)*a(3)+a(8)*a(6)];
return