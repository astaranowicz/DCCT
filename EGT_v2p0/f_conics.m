%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox  (EGT)  %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  f_conics    Plot of conics.
%
%  f_conics(A,color,passo,ampiezza)
%  A = a symmetric 3by3 matrix which defines the conic. A conic, in fact,
%      is understood as the locus of 2D points [x,y]' which satisfy the
%      following quadratic equation:
%             a*x^2 + b*x*y + c*x + d*y^2 + e*y + f = 0
%      This quadratic form can be rewritten in the matrix form U'AU = 0
%      where A is a symmetric 3by3 matrix 
%                           [  a    b/2   c/2 
%                             b/2    d    e/2
%                             c/2   e/2    f  ]
%      and U = [x,y,1]'. 
%  fig = integer to define the figure into which drawing.
%  color = character string which defines the colour used to draw.
%
%  Descr: 
%  -----  This function plots the conic which is defined by the matrix A
%         (if the conic is not singular).
%
%  f_conics(A,fig,color,passo,ampiezza) plots the conic (if no singular) defined in A.
%  flag = f_conics(A,fig,color) returns in flag an integer,
%         flag = 0  -->  regular conic
%         flag = 1  -->  singular conic: two imaginary lines with a 
%                                        real intersection
%         flag = 2  -->  singular conic: two transverses
%         flag = 3  -->  singular conic: two parallel lines
%

function flag = f_conics(A,color,passo,ampiezza)


if nargin==2,
 passo=1/100;
 ampiezza=10;
elseif nargin==3,
  ampiezza=10;
end;

if ((A(2,2)-A(1,1))==0)|((abs(A(2,2)-A(1,1)))<2^-52)
  phi = pi/4;
else
  phi = (atan(2*A(1,2)/(A(2,2)-A(1,1))))/2;
end;  
rotphi = rotoz(-phi);
% A1 = rotphi*A*rotphi';
A1 = rotphi'*A*rotphi;

detA = det(A);
delta = det(A(1:2,1:2));

% R01 = inv(rotphi);
R01 = rotphi;

% se il determinante ha un valore al di sotto di una certa soglia lo metto a 0
%   if abs(detA)<=2^-52%52
%       detA=0;
%   end
%   if abs(delta)<=2^-52%52
%       delta=0;
%   end
  
if detA~=0
    flag = 0;
    if delta~=0
        u01 = -(A1(1,3)/A1(1,1));
        v01 = -(A1(2,3)/A1(2,2));
    else
        if (A1(2,2)==0)|(abs(A1(2,2))<2^-52) 
            u01 = -(A1(1,3)/A1(1,1));
            v01 = (-A1(1,1)*u01^2-2*A1(1,3)*u01-A1(3,3))/(2*A1(2,3));
        elseif (A1(1,1)==0)|(abs(A1(1,1))<2^-52)
            v01 = -(A1(2,3)/A1(2,2));
            u01 = (-A1(2,2)*v01^2-2*A1(2,3)*v01-A1(3,3))/(2*A1(1,3));
        elseif A1(1,1)>=A1(2,2)
            u01 = -(A1(1,3)/A1(1,1));
            v01 = (-A1(1,1)*u01^2-2*A1(1,3)*u01-A1(3,3))/(2*A1(2,3));
        else
            v01 = -(A1(2,3)/A1(2,2));
            u01 = (-A1(2,2)*v01^2-2*A1(2,3)*v01-A1(3,3))/(2*A1(1,3))
        end
    end
 % le condizioni di maggioranza o minoranza servono per evitare che si 
 % inchiodi laddove dovesse trovare A1(1,1) e A1(2,2) entrambi diversi da 0
 % (può succedere?)
    
    A2 = [    A1(1,1)                 0                  A1(1,1)*u01+A1(1,3);
                 0                 A1(2,2)               A1(2,2)*v01+A1(2,3);
          A1(1,1)*u01+A1(1,3)  A1(2,2)*v01+A1(2,3)  A1(1,1)*u01^2+A1(2,2)*v01^2+A1(3,3)+2*A1(1,3)*u01+2*A1(2,3)*v01];
    if (delta>0) %ELLISSE
        flag = -1;
        [u2,v2] = f_ellipse(A2,passo,ampiezza);
        u1 = u2+u01;
        v1 = v2+v01;
        u = R01(1,:)*[u1;v1;ones(1,length(u1))];
        v = R01(2,:)*[u1;v1;ones(1,length(u1))];
        plot(u,v,color);

    elseif (delta<0) %HYPERBOLA
        [u2,v2] = f_hyperbola(A2,passo,ampiezza);
        u1 = u2+u01;
        v1 = v2+v01;  
        u = R01(1,:)*[u1;v1;ones(1,length(u1))];
        v = R01(2,:)*[u1;v1;ones(1,length(u1))];
          plot(u,v,color);
        flag = -1;
        u1neg = -u2+u01; 
        uneg = R01(1,:)*[u1neg;v1;ones(1,length(u1))]; 
        vneg = R01(2,:)*[u1neg;v1;ones(1,length(u1))];
        plot(uneg,vneg,color);

    else  %PARABOLA
        [u2,v2] = f_parabola(A2,passo,ampiezza);
        u1 = u2+u01;
        v1 = v2+v01;
        u = R01(1,:)*[u1;v1;ones(1,length(u1))];
        v = R01(2,:)*[u1;v1;ones(1,length(u1))];
          plot(u,v,color);         
        flag = -1;
    end
else
        disp('  EGT Warning: singular conic defined by the matrix')
end
    