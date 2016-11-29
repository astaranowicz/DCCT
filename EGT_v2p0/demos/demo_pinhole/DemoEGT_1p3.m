%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox  (EGT) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DemoEGT of Epipolar Geometry Toolbox
% ------------------------------------
%
% When run this demo places points in 3D space, places cameras and computes 
% the epipolar geometry. Epipolar lines are computed. Moreover the
% Fundamental matrix is also estimated with the EGT v1.3 algorithms and
% comparison on epipole estimation are reported.
%
% Last Update _ December 05
% by Gian Luca Mariottini
%

  close all; clear all;
  format long e;
    
% %Scene Points: Rectangular panels  
% %         NOTE = if the panels are added, due to the planarity of the scene the epipolar geometry estimation
% %         , as well known in literature, becomes singular and the Fundamental Matrix F is not 
% %         well estimated. Another feature point (last column above) has been be added to overcome 
% %         this singularity.
%X=0.1*[-5, 5,  5, -5,-5, 5, 5 , -5, 10, -5-1, 5-1,  5-10, -5-10, -5-10, 5-10, 5-10 , -5-10];
%Z=0.1*[ 5, 5, 15, 15, 5, 5,15 , 15, 10,   5-10, 5-10, 15-10, 15-10, 5-10, 5-10,15-10 , 15-10];
%Y=0.1*[15,15, 15, 15,25,25, 25, 25, 30,  15+10,15+10, 15+10, 15+10,25+10,25+10, 25+10, 25+10];
%P=[X;Y;Z];
  


%Scene Points: Random Points  
P=f_3Drandpoint(20,[0 2 .5],.6);

%%% PANEL
%%%
%Xpan=10^-1*[-5 5 5 -5];
%Ypan=[ 0 0 0  0];
%Zpan=10^-1*[-5 -5 5 5];
% 
% ang1=0;
% for i=1:length(Xpan),
%     Pan1([1:3],i)=rotoz(ang1)*[Xpan(i);Ypan(i);Zpan(i)];
% end;
% Xp1=Pan1(1,:);
% Yp1=Pan1(2,:);
% Zp1=Pan1(3,:);
% 
% ang2=pi/6;
% for i=1:length(Xpan),
%     Pan2([1:3],i)=rotoz(ang2)*[Xpan(i);Ypan(i);Zpan(i)];
% end;
% Xp2=Pan2(1,:);
% Yp2=Pan2(2,:);
% Zp2=Pan2(3,:);
% 
% ang3=-pi/6;
% for i=1:length(Xpan),
%     Pan3([1:3],i)=rotoy(-ang3)*rotoz(ang3)*[Xpan(i);Ypan(i);Zpan(i)];
% end;
% Xp3=Pan3(1,:);
% Yp3=Pan3(2,:);
% Zp3=Pan3(3,:);
% 
% n1=length(Xp1);
% P1=[Xp1+0*ones(1,n1) ; Yp1+4.5*ones(1,n1) ; Zp1+1.5*ones(1,n1)]; %350 su Y   
% n2=length(Xp2);
% P2=[Xp2-1.5*ones(1,n2) ; Yp2+4.0*ones(1,n2)  ; Zp2+.5*ones(1,n2)];  
% n3=length(Xp3);
% P3=[Xp3+1.5*ones(1,n3) ; Yp3+4.0*ones(1,n3)  ; Zp3+1.0*ones(1,n3)];     
% P=[P1 P2 P3];

figure(1); hold on; axis equal; grid on

%% WORLD REFERENCE FRAME
  f_3Dwf('k',1.5); % World reference frame
  f_3Dwfenum(P,'k',0.03);
%% SCENE POINTS
  f_scenepnt(P,'r*',1); %EGT->Plot of scene pnts.
  
%% FIRST (Desired) CAMERA
  % Camera plot
  Kd=[700 0 640; 
      0 700 480; 
      0    0   1];
  Rd=rotoy(-pi/2)*rotox(0)*rotoz(0); 
  td=[2.0,1,0]';
  Hd=f_Rt2H(Rd,td);
  f_3Dframe(Hd,'g',0.6,'_{d}');
  f_3Dcamera(Hd,'g',0.15,2); 
  U_d=f_perspproj(P,Hd,Kd,1);
  ud=U_d(1,:);
  vd=U_d(2,:);
  % Imaging model
  figure(2); hold on;  title('Image Plane - Desired Camera (green)'); grid on;
  xlabel('u [pixels]');ylabel('v [pixels]');
  plot(ud,vd,'gO');
  
%% SECOND (Actual) CAMERA
  % Camera plot
  figure(1)
  Ka=[700 0 640; 0 700 480; 0    0   1];
  Ra=rotoy(0)*rotox(0)*rotoz(0);
  ta=[0,0,0]';
  Ha=f_Rt2H(Ra,ta);
  f_3Dframe(Ha,'r',0.6,'_{c}'); 
  f_3Dcamera(Ha,'r',0.15,2);
  U_a=f_perspproj(P,Ha,Ka,1); % ua is a 2xN matrix, with N = "columns of P"
  ua=U_a(1,:);
  va=U_a(2,:);
  view(45,34);
  % Imaging model
  figure(3); hold on; title('Image Plane - Current Camera (red)'); grid on;
  xlabel('u [pixels]');ylabel('v [pixels]');
  plot(ua,va,'rO');

% Epipoles and FUndamental Matrix Computation 
  [ea,ed,F]=f_epipole(Ha,Hd,Ka,Kd);
  figure(2); plot(ed(1),ed(2),'kX'); text(ed(1),ed(2),'Epipole');
  figure(3); plot(ea(1),ea(2),'kX'); text(ea(1),ea(2),'Epipole');
  
  sigma=.3;
  Ua=[ua+sigma*randn(1,length(ua(1,:)));
      va+sigma*randn(1,length(ua(1,:)))];
  Ud=[ud+sigma*randn(1,length(ua(1,:)));
      vd+sigma*randn(1,length(ua(1,:)))];

  %Fstim4=f_Festim(Ua,Ud,2); % Robust M-estimator by Torr
  %[la,ld]=f_epipline(Ua,Ud,Fstim4,1,3,2,'r');
  %figure(4);
  display('[EGT message]: Ready to compare estimation methods over 30 iteration..')

for j=1:30,
  Ua=[ua+sigma*randn(1,length(ua(1,:)));
      va+sigma*randn(1,length(ua(1,:)))];
  Ud=[ud+sigma*randn(1,length(ua(1,:)));
      vd+sigma*randn(1,length(ua(1,:)))];
%% Epipolar Geometry Estimation with EGT
tic
Fstim1=f_Festim(Ua,Ud,4,-5,0.01); % Linear estimation with no det(F)=0
t1(j)=toc;
  tic;
  Fstim2=f_Festim(Ua,Ud,4,-16,0.01); % Normalized 8-point algorithm
  t2(j)=toc;

tic;
Fstim3=f_Festim(Ua,Ud,5); % Minimization of Sampson Distance
t3(j)=toc;
tic;
Fstim4=f_Festim(Ua,Ud,6); % Robust M-estimator by Torr
t4(j)=toc;
%Epipole in Alg.1
  if Fstim1~=zeros(3,3),
      e1a=null(Fstim1);
      e1d=null(Fstim1');
      e1a=e1a/e1a(3);
      e1d=e1d/e1d(3);
  else
      display('EGT warning: EGT can not find epipoles with METHOD 1 (linear unconstrained alg.). An infty value is substituted');
      e1a=[inf inf 1]'; e1d=[inf inf 1]';
  end    
% Epipole in Alg.2
  if Fstim2~=zeros(3,3),
      e2a=null(Fstim2);
      e2d=null(Fstim2');
      e2a=e2a/e2a(3); e2d=e2d/e2d(3);
  else
      display('EGT warning: EGT can not find epipoles. infty value is substituted');
      e2a=[inf inf 1]'; e2d=[inf inf 1]';
  end    
% Epipole in Alg.3
  if Fstim3~=zeros(3,3),
      e3a=null(Fstim3);
      e3d=null(Fstim3');
      e3a=e3a/e3a(3); e3d=e3d/e3d(3);
  else
      display('EGT warning: EGT can not find epipoles. infty value is substituted');
      e3a=[inf inf 1]';
      e3d=[inf inf 1]';
  end    
% Epipole in Alg.4
  if Fstim4~=zeros(3,3),
      e4a=null(Fstim4);
      e4d=null(Fstim4');
      e4a=e4a/e4a(3); e4d=e4d/e4d(3);
  else
      display('EGT warning: EGT can not find epipoles. infty value is substituted');
      e4a=[inf inf 1]';
      e4d=[inf inf 1]';
  end    
% Epipole errors
err1a(j)=norm(e1a-ea);err1d(j)=norm(e1d-ed);
err2a(j)=norm(e2a-ea);err2d(j)=norm(e2d-ed);
err3a(j)=norm(e3a-ea);err3d(j)=norm(e3d-ed);
err4a(j)=norm(e4a-ea);err4d(j)=norm(e4d-ed);
end
% Plot of estimation results
figure(4); hold on;title('Epipole estimation error: comparison over N=30 iterations');
xlabel('Algs:  1-Lin.(Norm.8pnt) , 2-Iter.(Samp.) , 3-M-estim. (Torr)'); ylabel('Error Mean/Variance [pixels]');
ERROR_ESTIM=[mean(err2a) mean(err2d);
             mean(err3a) mean(err3d);
             mean(err4a) mean(err4d)]; 
bar([1:3] , ERROR_ESTIM);
legend('Error on epipole in <c> (pix)','Error on epipole in <d> (pix)')
errorbar(1-0.15, mean(err2a), std(err2a));
errorbar(1+0.15, mean(err2d), std(err2d));
errorbar(2-0.15, mean(err3a), std(err3a));
errorbar(2+0.15, mean(err3d), std(err3d));
errorbar(3-0.15, mean(err4a), std(err4a));
errorbar(3+0.15, mean(err4d), std(err4d));


figure(5); hold on;
bar([1:3], [mean(t2), mean(t3), mean(t4)]);
errorbar(1, mean(t2), std(t2));
errorbar(2, mean(t3), std(t3));
errorbar(3, mean(t4), std(t4));
xlabel('Algs:  1-Lin.(Norm.8pnt) , 2-Iter.(Samp.) , 3-M-estim. (Torr)'); ylabel('Mean/Variance [secs.]');


