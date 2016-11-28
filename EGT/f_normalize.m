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
    for i=1:npunti,
	    Usc(:,i)=U(:,i)-[cu,cv]';
    end,
    %Adesso che hai scalato i punti nel centroide calcolane la norma
    for i=1:npunti,
        NORMA(i)=norm(Usc(1:2,i));   
    end,
    NORMAtot=sum(NORMA)/npunti; 
    t=sqrt(2)/NORMAtot;
    T=[t,0,-t*cu;
       0,t,-t*cv;
       0,0,   1];   
    for i=1:npunti,
        Uo([1:3],i)=[U(1,i) ; U(2,i) ; 1];
    end;
    
    Unormsco=T*Uo;   %la T lavora sui punti non traslati nel centroide  
    Unormsc([1,2],:)=[Unormsco(1,:);Unormsco(2,:)];
    
    %% Decomment this block to see centroid and distance between trasformed
    %% point values
    %%
    %
    %
    %display('Mean distance from center is...')
    %for i=1:npunti,
    %    NORMAnorm(i)=norm(Unormsc(1:2,i));
    %end,
    %distanzammediaFINALE=sum(NORMAnorm)/npunti
    %centroideFINALE=[sum(Unormsc(1,:))/npunti;
    %                 sum(Unormsc(2,:))/npunti]
    
    