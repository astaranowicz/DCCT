%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       %%
%%  Epipolar Geometry Toolbox v1.1 (EGT) %%
%%                                       %%
%%%%%%%% DII- University of Siena %%%%%%%%%
%
%  [Unormsc,T]=f_normalize(U);
%
% Descr.=
% ------ This function performs a normalization of input points so that 
%        the centroid of the reference points s in (0,0)' and the RMS 
%        distance of the transformed points from the origin is equal to 
%        sqrt(2).
%
% Syntax=
% ------ U = "matrix of points in 3d scene" (see manual) (can be also in homogeneous notation)
%        T = "normalization matrix"
%        Unormsc = "normalized and scaled matrix of points"
%
% Gian Luca Mariottini - October/November 2003
%
function [Unormsc,T]=f_normalize(U);
%given a set of points U the normalization is performed
npunti=length(U(1,:));
%1] Centroid computation
    cu=sum(U(1,:))/npunti;
    cv=sum(U(2,:))/npunti;
%2] Scaling factor
    
	Usc([1:2],:)=U([1:2],:)-[cu;cv]*ones(1,npunti);
    
    %Adesso che hai scalato i punti nel centroide calcolane la norma    
    NORMAtot=sqrt(Usc(1,:)*Usc(1,:)' + Usc(2,:)*Usc(2,:)')/npunti; 
    
    t=sqrt(2)/NORMAtot;
    T=[t,0,-t*cu;
       0,t,-t*cv;
       0,0,   1];   
    
        Uo([1:3],:)=[U(1,:) ; U(2,:) ; ones(1,npunti)];
    
    
    Unormsco=T*Uo;   %la T lavora sui punti non traslati nel centroide  
    Unormsc([1,2],:)=[Unormsco(1,:);Unormsco(2,:)];
    
    