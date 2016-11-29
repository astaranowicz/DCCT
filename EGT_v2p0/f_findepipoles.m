
%%%%%%%%%%%%%%%%%%%%%%%
%%%%  subfunction  %%%%
%%%%%%%%%%%%%%%%%%%%%%%

function [qe1, qe2, e1, e2] = f_findepipoles(F,a,b,K);
% calcola gli epipoli o meglio la loro proiezione sulo specchio (sdr mirror) 
% F = coordinate del fuoco dell'altro specchio espresse nel sdr mirror 
% relativo allo specchio per il quale sto calcolando gli epipoli
% f = [e1,e2,qe1,qe2] : la funzione restituisce sia gli epipoli sullo
% specchio che le rispettive proiezioni sul piano immagine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K = [K(1,1) K(1,3) K(1,2); %
%      K(2,1) K(2,3) K(2,2); %
%      K(3,1) K(3,3) K(3,2)];%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Rcam = eye(3);
  tcam = [0 0 2*sqrt(a^2+b^2)]';


F1 = [0,0,0]'; F2 = F;
xa = F1(1); ya = F1(2); za = F1(3); 
xb = F2(1); yb = F2(2); zb = F2(3); 
l = xb-xa ; m = yb-ya ; n = zb-za ;
e = sqrt(a^2+b^2) ;
% Ax^2 + B*x + C = 0;
A = (b^2*n^2/l^2)-a^2*(1+m^2/l^2) ;
B = 2*a^2*(xa*m^2/l^2-za*m/l)+2*b^2*((ya+e)*n/l-xa*n^2/l^2) ;
C = b^2*(xa^2*n^2/l^2+(ya+e)^2-2*xa*(ya+e)*n/l)+a^2*(2*za*xa*m/l-xa^2*m^2/l^2-za^2-b^2) ; 
x1 = (-B-sqrt(B^2-4*A*C))/(2*A) ;
x2 = (-B+sqrt(B^2-4*A*C))/(2*A) ;
y1 = (x1-xa)*m/l+ya ;
y2 = (x2-xa)*m/l+ya ;
z1 = (y1-ya)*n/m+za ;
z2 = (y2-ya)*n/m+za ;

if F(1)~=0,
    e1 = [x1,y1,z1]';
    e2 = [x2,y2,z2]';
% e1 ed e2 sono le coordinate nel sdr mirror degli epipoli trovati come
% intersezione tra la retta congiungente i due fuochi degli specchi e lo
% specchio al quale si riferisce il sdr mirror usato
    e1cam = Rcam*e1+tcam; % e1cam = coordinate di e1 nel sdr camera
    qe1 = K*(1/e1cam(3))*e1cam;
    e2cam = Rcam*e2+tcam;
    qe2 = K*(1/e2cam(3))*e2cam;
else %Planar motion case...
    %
    % Remember that the conic equation is (z+e)^2/a^2 - (x^2+y^2)/b^2 =1. 
    %
    e1 = [0, sqrt(b^2*(e^2/a^2-1)),0]';
    e2 = [0,-sqrt(b^2*(e^2/a^2-1)),0]';
% e1 ed e2 sono le coordinate nel sdr mirror degli epipoli trovati come
% intersezione tra la retta congiungente i due fuochi degli specchi e lo
% specchio al quale si riferisce il sdr mirror usato
    e1cam = Rcam*e1+tcam; % e1cam = coordinate di e1 nel sdr camera
    qe1 = K*(1/e1cam(3))*e1cam;
    e2cam = Rcam*e2+tcam;
    qe2 = K*(1/e2cam(3))*e2cam;
end    