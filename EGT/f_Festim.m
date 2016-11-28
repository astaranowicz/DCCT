%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.1 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% function f_Festim(Ua,Ud,algorithm)
%
% Descr: 
% -----  This function performs the estimation of fundamental matrix
%        F from a set of given correspondent points Ua and Ud. 
%        Several algorithms are here implemented.
%        (see manual for details)
%
% Syntax:
% ------  Ua=matrix 2xn , where n is the number of feature points (see the
%         manual).
%         "algorithm"= 1,  "use linear estimation with no det(F)=0 constraint"
%                    = 2,  "Normalized 8-point algorithm"
%                    = 3,  "Minimization of Sampson Distance"
%                    = 4,  "Robust Torr's M-estimator with outliers
%                    removal"
% Author:
%     Gian Luca Mariottini 
% Last update:
%     December 2005
function F=f_Festim(Ua,Ud,algorithm);

if algorithm==1,
    for i=1:length(Ua(1,:)),
        A(i,[1:9])=[Ud(1,i)*Ua(1,i) , Ud(1,i)*Ua(2,i) , Ud(1,i) , Ud(2,i)*Ua(1,i) , Ud(2,i)*Ua(2,i) , Ud(2,i) , Ua(1,i) , Ua(2,i) , 1];
    end;
     if rank(A)<8,
        display('EGT error: Fundamental Matrix estimation "error" ! =>  rank(A)<8  <=');
     elseif rank(A)==8,
        %display('EGT communication: Fundamental Matrix "linear estimation" in progress...');
        f=null(A);
     elseif rank(A)==9,
        display('EGT communication: "8 point algorithm" for F estimation in progress...');
        %
        %==> 8 point algorithm 
        %
        %The least sqaures solution for f is the singular vector corresponding to the smallest  
        %singular value of A (i.e the last column of V in SVD(A)=UDV'
        %D contains, in descending order, the singualr values of matrix A;
        [U,D,V]=svd(A);
        f=V(:,length(V(1,:)));
     end;
     if rank(A)>=8,
        sf=5;
        f=f_roundn(f,-sf); %Arrotondo alla settima cifra decimale
        F=[f(1) f(2) f(3);
           f(4) f(5) f(6);
           f(7) f(8) f(9)];       
     else
        display('EGT error: f_Festim has too few points to estimate the Epipolar Geometry!')
        F=zeros(3,3);
        return
    end;   
     % Rank Check
     if rank(F)~=2,
        display('EGT Warning: EGT can not find the Fundam.Matrix with Linear Unconstrained Alg.!');
        display('             the matrix F=0 is given in output.')    
        F=zeros(3,3);
     end;    
elseif algorithm==2,
     %1step] Normalization of feature points
        [Uanorm,Ta]=f_normalize(Ua);
        [Udnorm,Td]=f_normalize(Ud);
    for i=1:length(Ua(1,:)),
        A(i,[1:9])=[Udnorm(1,i)*Uanorm(1,i) , Udnorm(1,i)*Uanorm(2,i) , Udnorm(1,i) , Udnorm(2,i)*Uanorm(1,i) , Udnorm(2,i)*Uanorm(2,i) , Udnorm(2,i) , Uanorm(1,i) , Uanorm(2,i) , 1];
    end;
     % ----------------------------
     % Normalized 8-point algorithm    (in ["Multiple View Geometry"
     % ----------------------------     R.Hartley,A.Zisserman, 2000-pg.265-266])  
     %      
      if rank(A)<7,
        display('EGT error: Fundamental Matrix estimation "error" ! =>  rank(A)<8  <=');
      else
        % 1] Linear solution
           [U,D,V]=svd(A);  %Singular decomposition of A (in D all eigenvalues of A in descending order are reported)
           f=V(:,length(V(1,:))); %vector corresponding to the samllest eigenvalue
           sf=5;
           f=f_roundn(f,-sf); %Round to the "sf" decimal digit
           Func=[f(1) f(2) f(3);
                 f(4) f(5) f(6);
                 f(7) f(8) f(9)];%Unconstrained matrix
        % 2] Denormalization
           Funcden=Td'*Func*Ta;
        % 3] Constraint enforcement
           [Uf,Df,Vf]=svd(Funcden);
           rf=Df(1,1);
           sf=Df(2,2);
           Ff=Uf*diag([rf,sf,0])*Vf';
           F=Ff;        
       end;
     % Rank Check
     if rank(F)~=2,
        display('EGT Warning: EGT can not find the Fundam.Matrix with 8-point Algorithm!');
        display('             the matrix F=0 is given in output.');        
        F=zeros(3,3);
     end;     
elseif algorithm==3,
     F=f_grad(Ua,Ud); %include normalization
 elseif algorithm==4,
      [Uanormsc,Ta]=f_normalize(Ua); %Actual points
      [Udnormsc,Td]=f_normalize(Ud); %Desired points
      Mn=[Udnormsc; Uanormsc];
      %Mn=[Ud; Ua];
      % F estimation
      errorestimaF=10^-1;
      iterazioni=10;
      [F,ga]=f_FestimM_byTorr(Mn,iterazioni,errorestimaF);       
      F=Td'*F*Ta; %(d',a) De-normalize fundamental matrix    
 end     
 
 
function [F,g]=f_FestimM_byTorr(M,n,error);
   % Check for errors
   if (size(M,1)~=4) | (size(M,2)<8),
      disp('EGT Error: error in parameter declaration')
   else
   Fant=[0 0 0; 0 0 0; 0 0 0];
   U=[];
   for i=1:size(M,2),
  	   U=[U; M(1,i)*M(3,i) M(1,i)*M(4,i) M(1,i) M(2,i)*M(3,i) M(2,i)*M(4,i) M(2,i) M(3,i) M(4,i) 1];
   end
   g=ones(1,size(M,2));
   w=ones(1,size(M,2));
   for iter=1:n,      
   w2=w.*g;
     for i=1:size(M,2),
       U2(i,:)=U(i,:).*w2(i);         
     end       
   [V,D]=eig(U2'*U2);      
   [min_val,ind_i]=min(sum(D));
   f=V(:,ind_i);	
   F=f_unstacked(f,3)';      
   F=f_forcerank2(F);
   if max(max(abs(F-Fant)))<error,
     break;
   end
   Fant=F;      
   for i=1:size(M,2),
         tr(i)=[M(1:2,i) ; 1]'*F*[M(3:4,i) ; 1];
         x1=M(1,i);
         y1=M(2,i);
         x2=M(3,i);
         y2=M(4,i);
         rx2=F(1,1)*x1+F(2,1)*y1+F(3,1);
         ry2=F(1,2)*x1+F(2,2)*y1+F(3,2);
         rx1=F(1,1)*x2+F(1,2)*y2+F(1,3);
         ry1=F(2,1)*x2+F(2,2)*y2+F(2,3);
         w(i)=sqrt(1/(rx2^2+ry2^2+rx1^2+ry1^2));
         dt(i)=w(i)*tr(i);
      end;
      s_i=median(abs(dt))/0.6745;
      for i=1:size(M,2),
         if abs(dt(i))<3*s_i, g(i)=s_i/abs(dt(i));
         elseif abs(dt(i))<s_i, g(i)=1;
         else  g(i)=0; end
      end
   end
end

function [M]=f_forcerank2(Mi);
if (size(Mi,1)~=3) | (size(Mi,2)~=3),
   disp('EGT Error: incorrect numebr of parameters')
else
   [U,S,V] = svd(Mi);
   S(3,3)=0;
   M = U*S*V';
end