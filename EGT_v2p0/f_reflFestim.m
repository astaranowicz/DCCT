%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.1 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
% 
% function f_Festim(U,U',algorithm)
%
% Descr: 
% -----  This function performs the estimation of fundamental matrix
%        F from a set of given correspondent points Ua and Ud. 
%        Several algorithms are here implemented.
%        (see manual for details)
%
% Syntax:
% ------  U=matrix 2xn , where n is the number of feature points (see the manual).
%         "algorithm"= 1,  "use linear estimation with no det(F)=0 constraint"
%                    = 2,  "Normalized 8-point algorithm"
%                    = 3,  "Minimization of Sampson Distance"
%                    = 4,  "Robust Torr's M-estimator with outliers removal"
% Author:
%     Gian Luca Mariottini 
% Last update:
%     March 2008
function F=f_Festim(U,Up,algorithm,roundf);
global ga;
ga=zeros(1,length(U(1,:)));
if nargin==3,
    roundf=-16;
end

if algorithm==1,
    n=length(U(1,:));
    % Create the matrix A 
    A=zeros(n,9);
    for i=1:n,
        A(i,[1:9])=[Up(1,i)*U(1,i) , Up(1,i)*U(2,i) , Up(1,i) , Up(2,i)*U(1,i) , Up(2,i)*U(2,i) , Up(2,i) , U(1,i) , U(2,i) , 1];
    end;
     if rank(A)<8,
        display('EGT error: Fundamental Matrix estimation "error" ! =>  rank(A)<8  <=');
        display('           f_Festim has too few points to estimate the Epipolar Geometry!')
        F=zeros(3,3);
        return
     elseif (rank(A)==8)||(rank(A)==9),
        display('EGT communication: "8 point algorithm" for F estimation in progress...');
        %
        %==> 8 point algorithm 
        %
        %The least sqaures solution for f is the singular vector corresponding to the smallest  
        %singular value of A'*A (i.e the last column of V in SVD(A'*A)=UDV'
        %D contains, in descending order, the singualr values of matrix A'*A (Ziss.pag 90);
        [U,D,V]=svd(A'*A);
        f=V(:,length(V(1,:)));
        % Rounding and imposing the det(F)=0 contraint
        f=f_roundn(f,roundf); 
        F=[f(1) f(2) f(3);
           f(4) f(5) f(6);
           f(7) f(8) f(9)];   
        %F=f_forcerank2(F);  
    end; 
elseif algorithm==2,
     % ----------------------------
     % Normalized 8-point algorithm  
     % ----------------------------     
     % (pag.282 - Zisserman -2nd. ed.
     
     % 1st step] Normalization of feature points     
        [U_n ,T ] = f_normalize(U);
        [Up_n,Tp] = f_normalize(Up);
        
     % 2nd step] Obtain the matrix F   
        Fhat_p=f_Festim(U_n,Up_n,1);    

     % 3rd step] Denormalization
        F = Tp'*Fhat_p*T;    
        n=norm(F,'fro');
        if n~=0,
              F=F./n;
        end
%% 3) ITERATIVE Minimization with Sampson distance + Fund. matrix. param.    
elseif algorithm==3,
    % 1- Linear solution
     F0=f_Festim(U,Up,2); % [Up([1:2],1);1]'*F0*[U([1:2],1);1] %OK!    
    % 2- Nonlinear refinement 
      options1  = optimset('Algorithm','levenberg-marquardt','TolX',1e-6,'TolFun',1e-11,'Display','off'); %6 e 11
      global U Up;
      f0= f_stack(F0);
      a=fminunc(@sampson,f0,[],[],options1);
%     % 3- Solution 
%      F=unstkc(a,3,3);
%Mn=[U; Up];
%[f,f_sq_errors] = torr_estimateF(Mn(:,:)', 1, [], 'non_linear',1,f0);
            
          F0=unstkc(a,3,3);
          F=F';
          
     n=norm(F,'fro');
        if n~=0,
              F=F./n;
        end
      
%% 4) Torr Algorithm      
elseif algorithm==4,
      [U_normsc,T]=f_normalize(U); %Actual points
      [Up_normsc,Tp]=f_normalize(Up); %Desired points
      Mn=[U_normsc; Up_normsc];
      % Parameters
      error_estim_F_Torr=10^-5;
      iterations_Torr=500;
      [F,ga]=f_FestimM_byTorr(Mn,iterations_Torr,error_estim_F_Torr,Tp,T);
%         Thresh_ga_outli=(min(ga)+median(ga))/2;
%         ind_inliers=find(abs(ga)>=Thresh_ga_outli);
%         if length(ind_inliers)<7,
%           Thresh_ga_outli=(min(ga)+median(ga))/2;
%           ind_inliers=find(abs(ga)>=Thresh_ga_outli);
%         end 
%         [f,f_sq_errors] = torr_estimateF(Mn(:,ind_inliers)', 1, [], 'non_linear',1,f_stack(F)); 
%         F=unstkc(f,3,3)
      F=F';
      F=Tp'*F*T; %(d',a) De-normalize fundamental matrix   
      n=norm(F,'fro');
        if n~=0,
              F=F./n;
        end
             
%% 5) LMedS Algorithm using Eigenvalues    
elseif algorithm==5, 
      %% Hartley's point normalization
        [U_normsc,T]=f_normalize(U); %Actual points
        [Up_normsc,Tp]=f_normalize(Up); %Desired points
        Mn=[U_normsc;
            Up_normsc];
      %% F estimation
        b=8;%7; %bucket
        p=0.9; %Probability that there is an F without outliers
        r=0.25; %Outlier_ratio
        [F,ga]=funmatLMedSeig(Mn,b,p,r); 
      %% Outliers removal              
          Thresh_ga_outli=(min(ga)+median(ga))/2;
          ind_inliers=find(abs(ga)>=Thresh_ga_outli);
          if length(ind_inliers)<=7,
              ind_inliers=find(abs(ga)>0);
          end
      %% Refinement using all inliers
        [F,ga]=funmatLMedSeig(Mn(:,ind_inliers),b,p,r);   
          
      F=F';
      F=Tp'*F*T; %(a',d) De-normalize fundamental matrix
              
      
      F=f_forcerank2(F);
      n=norm(F,'fro');
        if n~=0,
              F=F./n;
        end
              
elseif algorithm==6,
      %% Hartley's point normalization
      [U_normsc,T]=f_normalize(U); % Actual points
      [Up_normsc,Tp]=f_normalize(Up); % Desired points
      Mn=[U_normsc;
          Up_normsc];
      %% Robust MAPSAC F-estimator
        %-1- Initial estimate
            global inlier_index;
            [f, f_sq_errors, n_inliers,inlier_index,Fmapsac] = torr_estimateF( Mn', 1, [round(log(1-0.99)/log(1-(1-0.25)^7)) 0.0039], 'mapsac', 0);
        %-2- Nonlinear minimization refinement
            [f,f_sq_errors] = torr_estimateF(Mn(:,inlier_index)', 1, [], 'non_linear',1,f);
            
          F0=unstkc(f,3,3);
          F0=F0';  
          F=Tp'*F0*T;%(a',d) De-normalize fundamental matrix   
          
        %-3- Adaptation           
          F=f_forcerank2(F);      
end
%
 
 
function [F,g]=f_FestimM_byTorr(M,n,error,Td,Ta);
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
   
   % De-normalization
   %F=Td'*F*Ta;
   % Force rank
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
   disp('EGT Error: incorrect number of parameters')
else
   [U,S,V] = svd(Mi);
   S(3,3)=0;
   M = U*S*V';
end

% Distanza di Sampson (equazione 18)
function c=sampson(fundmat)
    global U Up nCols;
    fundmat;
    	Fe=unstkc(fundmat,3,3);
        Fe2=Fe;%f_forcerank2(Fe);
        f2=f_stack(Fe2);
    Fe2;
    %display('fine')
        if (rank(Fe2)==1),
            % In this case a parameterization must be used, as shown in 
            % "Parametrizations of the Essential and the Fundamental Matrix
            % based on Householder Transformations "
            % by Fabian Wenzel and Rolf-Rainer Grigat
            display('rank=1')
            % The matrix H must be built (pg.392 MaSKS)
            e2=null(Fe2);
            ind_inf = find (e2(3,:)==0);
            e = e2(:,ind_inf);
            e
            pause
            % Construction of matrix H
            alpha=e(1);
            beta=e(2);
            H=1/2*[1+beta^2-alpha^2, -2*alpha*beta    ,    -2*alpha;
                    -2*alpha*beta  , 1- beta^2+alpha^2,    -2*beta;
                        2*alpha    ,   2*beta         , 1-beta^2-alpha^2];
                    
            alpha1=alpha;
            beta1=beta;
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