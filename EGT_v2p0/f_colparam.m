%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox  (EGT) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function U=f_colparam(I,flag,type);
%
% where U=[mr mg mb vr vg vb sigma_r sigma_g sigma_b];
%
% Descr.
% -----
%   Computation of color parameters for color filtering.
%   The function take in input a image "I" and enable the user to click in a
%   (type=1) or to crop on a (type=2) image section. Pixel will be used to
%   compute mean and variance.
%
% Synt.
% ----
%   I= rgb image (acquired via imread of MATLAB)
%   flag = use adaptive or not technique for sigma computation
%   
%   type = 1 --> clicking
%   type = 2 --> croping
%   flag = 1 --> non adattativo
%   flag = 2 --> adattativo
%
% Author:
%     Gian Luca Mariottini 
% Last update:
%     June 2004
function U=f_colparam(I,flag,type)

sigma_r=2.5;
sigma_g=2.5;
sigma_b=2.5;
%  I=vfm('grab',1);
% 
% %Quale algoritmo per selezionare il colore ?
%  type=1;
%  flag=2;
% %%%
    
if flag==2,
    %Adaptative
    display('Select the area....');
    [I2,rect]=imcrop(I);
    rectangle('Position',rect,'EdgeColor','w');
    display('...select the color!');
end;

%
% BASIC COLOR FILTERING (flag==1)
%
figure(1)
imshow(I)
if type == 1
  [Yt,Xt]=ginput;
  for i =1:length(Xt),
    val=I(round(Xt(i)),round(Yt(i)),:);
    Int(:,i)=val;
  end;
  Int=double(Int);
  
  mr=mean(Int(1,:));
  mg=mean(Int(2,:));
  mb=mean(Int(3,:));
  vr=var(Int(1,:));
  vg=var(Int(2,:));
  vb=var(Int(3,:));
else
  I2 = imcrop(I);
  Int=double(I2);
  mr=mean(mean(Int(:,:,1)));
  mg=mean(mean(Int(:,:,2)));
  mb=mean(mean(Int(:,:,3)));
  vr=mean(var(Int(:,:,1)));
  vg=mean(var(Int(:,:,2)));
  vb=mean(var(Int(:,:,3)));
end;

  mv=[mr mg mb vr vg vb];

%
% ADAPTATIVE COLOR FILTERING
%
if flag==2,
%--> CICLO DI FILTRAGGIO <--
Tend=20,
Area_det(1)=rect(3)*rect(4);
for i=2:Tend;
     %I=vfm('grab',1);
     %figure(1);
     
     %1. Usa i parametri trovati solo cliccando sul colore
        I_filt=f_colfilt(I,mv(1),mv(4),mv(2),mv(5),mv(3),mv(6),sigma_r,sigma_g,sigma_b);
     %2. Per capire quale colore hai selezionato
        [massimo_media,indice_col_media]=max([mv(1),mv(2),mv(3)]);
        [minimo_varianza,indice_col_var]=min([mv(4),mv(5),mv(6)]);
        if (indice_col_media==1),
           colore_selez='red';
        elseif (indice_col_media==2),
           colore_selez='ver';
        elseif (indice_col_media==3),
           colore_selez='blu';
        end;

     %3. Vedi se la varianza del colore ci permette di filtrare roba buona o fuori dalla regione di effettivo interessse
        Ibw_fill = imfill(I_filt,[1 1]);
        Idouble=double(~Ibw_fill);
        J = medfilt2((Idouble),[5 5]);
        Jclose = J; %imfill(J,'hole'); %CHiude l'immagine
        C=f_regionprops(Jclose,'Centroid','BoundingBox','Area');
        if isempty(C),
           corr_x=corr_old;
           corr_y=corr_old;
        else
            Cx=C.Centroid(1);
            Cy=C.Centroid(2);
            BBx=C.BoundingBox(1);
            BBy=C.BoundingBox(2);
            BBw=C.BoundingBox(3);
            BBh=C.BoundingBox(4);
            Area=BBw*BBh;
            %Area=C.Area;
        end;
        %Elimina outliers
        if ((Area>Area_det(i-1)+2000)&(Area<Area_det(i-1)-2000)),
            Area_det(i)=Area_det(i-1);
        else
            Area_det(i)=Area;
        end;
        
        corr_old=C;
        corr_old=C;  
     
   %4. Adatta la varianza a seconda del colore  
      Diff2(i)=(Area_det(i)-rect(3)*rect(4));  
   %5. Quale funzionale adottare...e adottalo  
      if colore_selez=='ver',
          sigma_vecchio=sigma_g;
      elseif colore_selez=='red',
          sigma_vecchio=sigma_r;
      elseif colore_selez=='blu',
          sigma_vecchio=sigma_b;
      end,
  %Se si in un certo bound del box iniziale non cambiare nulla 
      soglia=250;
    if abs(Diff2(i)) > soglia,     
        %Dividi=2000;
        %sigma_mod(i)=sigma_vecchio-tanh((Diff2)/Dividi);
        funzionale=2;
        if funzionale==1,
            Dividi=2000;  
            Aggiorna(i)=-tanh((Diff2(i))/Dividi);
        elseif funzionale==2,
            Dividi=300000;  
            Aggiorna(i)=sign(Diff2(i))*10^-1*(1-exp(-Diff2(i)^2/Dividi));
            %pause
        end,
        %AGGIORNA LA VARIANZA
        sigma_mod(i)=sigma_vecchio-Aggiorna(i);
    else % Nel caso in cui sei DENTRO LA SOGLIA....non aggiornare sigma_mod
        sigma_mod(i)=sigma_vecchio;
        Aggiorna(i)=0;
    end;
    
    %6. Aggiorna i valori
        if colore_selez=='ver',
          sigma_g=sigma_mod(i);
          sigma_r=sigma_r-Aggiorna(i)/2;
          sigma_b=sigma_b-Aggiorna(i)/2;        
        elseif colore_selez=='red',
          sigma_r=sigma_mod(i);
          sigma_g=sigma_g-Aggiorna(i)/2;
          sigma_b=sigma_b-Aggiorna(i)/2;  
        elseif colore_selez=='blu',
          sigma_b=sigma_mod(i);
          sigma_r=sigma_r-Aggiorna(i)/2;
          sigma_g=sigma_g-Aggiorna(i)/2;   
        end,
  
  %figure(1)
  %imshow(J)
  %hold on
  %rectangle('Position',[BBx,BBy,BBw,BBh],'EdgeColor','g')
  %rectangle('Position',[BBx,BBy,rect(3),rect(4)],'EdgeColor','r')
  %pause(0.01)
  %hold off
end;

figure
plot(sigma_mod)
title('Behavior of adaptive variance w.r.t time')

figure
plot(Area_det,'g')  
hold on
plot(ones(1,Tend)*rect(3)*rect(4),'r--')
title('Area of selected image')
end;

U=[mv sigma_r sigma_g sigma_b];