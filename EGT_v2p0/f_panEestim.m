function [E,gE]=f_panEestim(M,n,error)
%Each column of M contains the mirror projections X_M of the 3-D point P in both actual and desired views.
%M=[X_Ma; 
%   X_Md]

%Controlla che siano dati in ingresso almeno le 6 coordinate nei 2 CCD
%ed almeno 8 punti nella scena
if (size(M,1)~=6) | (size(M,2)<8),
   disp('Error: parameters incorrect')
else
   v2=ones(1,size(M,2));
   gE=ones(1,size(M,2));
   
   Eant=[0 0 0; 0 0 0; 0 0 0];
   A=[];
 %Crea la matrice A tale che Ae=0
	for i=1:size(M,2),
        A=[A;(kron(M(1:3,i),M(4:6,i)))'];
	end
    
  %Inizio iterazioni
    for iter=1:n,      
      ww=v2.*gE;
      for i=1:size(M,2),
   	   AA(i,:)=A(i,:).*ww(i);         
      end      
      
      [V,D]=eig(AA'*AA);
         
      [minim,i]=min(sum(D));
	  e=V(:,i);
	
      E=[e(1) e(2) e(3); 
         e(4) e(5) e(6); 
         e(7) e(8) e(9)];
      
      %Forzo E ad avere rango 2
      E=f_forcerank2(E);
      
      if max(max(abs(E-Eant)))<error,
         break;
      end
      Eant=E;
      
       for i=1:size(M,2),
         r(i)=[M(1:3,i)]'*E*[M(4:6,i)];
         x1=M(1,i);
         y1=M(2,i);
         z1=M(3,i);
         x2=M(4,i);
         y2=M(5,i);
         z2=M(6,i);
         rx2=E(1,1)*x1+E(2,1)*y1+E(3,1)*z1;
         ry2=E(1,2)*x1+E(2,2)*y1+E(3,2)*z1;
         rz2=E(1,3)*x1+E(2,3)*y1+E(3,3)*z1;
         rx1=E(1,1)*x2+E(1,2)*y2+E(1,3)*z2;
         ry1=E(2,1)*x2+E(2,2)*y2+E(2,3)*z2;
         rz1=E(3,1)*x2+E(3,2)*y2+E(3,3)*z2;
         w(i)=sqrt(1/(rx2^2+ry2^2+rx1^2+ry1^2));
         d(i)=v2(i)*r(i);
      end;
      
      sig=median(abs(d))/0.6745;
      
      for i=1:size(M,2),
         if abs(d(i))<sig,
            gE(i)=1;
         elseif abs(d(i))<3*sig,
            gE(i)=sig/abs(d(i));
         else
            gE(i)=0;
         end
      end
   end %Fine iterazioni
end %Fine controllo parametri ingresso

function [M]=f_forcerank2(Mi);
if (size(Mi,1)~=3) | (size(Mi,2)~=3),
   disp('EGT Error: incorrect numebr of parameters')
else
   [U,S,V] = svd(Mi);
   S(3,3)=0;
   M = U*S*V';
end