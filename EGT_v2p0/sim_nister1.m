%% Initialize the environment
clear all;
close all;
clc;

figure(1);
hold on;
grid on;
title('3D setup');

axis equal;
view(-18,72);

%% Parameter
d_ch=1;  % Adjustable parameter for camera-head distance: if d_ch=0 then all the camera coincide with <h>

%% Create the 3D points
    % P=[P_p, Ptab] is a "n X 4" matrix
    % where P_p (nx3) contains the n 3D point coordinates
    % and P_tab (nx1) contains the integer representing the index f the
    %       corresponding image point.
      P_p = f_3Drandpoint(3,[0 0 1],1);  % n x 3 matrix
      P_tab = [1:length(P_p(1,:))];
      P = [P_p', P_tab'];
      f_scenepnt(P_p,'r*',1);
      f_3Dwfenum(P_p,'k',.05);
      
%% Create the Helmet <h>
      th_m = [0, -2, .5]';   % Pose of <h> in the world MATLAB frame <M>
      text(th_m(1)-.1,th_m(2),th_m(3),'<H>');
      Rh_m = eye(3); % This ROTATION is the rotation FROM <w> TO <h>
      Hh_m = [Rh_m   , th_m;
              [0 0 0], 1];
    % Pose of the 1st camera (wrt <h>)
      K(:,:,1) = [2000    0 640; 
                  0    2000 480; 
                  0       0   1];
      Rc_h(:,:,1) = rotoy(0)*rotox(0)*rotoz(0); 
      tc_h(:,1) = d_ch*[0, 0, .2]';
    % Pose of the 2nd camera (wrt <h>)   
      K(:,:,2) = K(:,:,1);
      Rc_h(:,:,2) = rotoy(-pi/6); 
      tc_h(:,2) = d_ch*[-.2, 0, .1]';
    % Pose of the 3rd camera (wrt <h>)
      K(:,:,3) = K(:,:,1);
      Rc_h(:,:,3) = rotoy(pi/6); 
      tc_h(:,3) = d_ch*[ .2, 0, .1]';     
    % ...add as many you want ;-)    
      
%% 3D Plot the helmet <h>   
      for i=1:3,
          Hc_h([1:4],[1:4],i)=[Rc_h(:,:,i), tc_h(:,i);
                               zeros(1,3) ,     1    ];
      end;
      %Plot the cameras on the helmet
          f_3Dhelmet(Hh_m, Hc_h, .1, 'r', 1);
      % Plot the helmet frame
          hold on
          Hh = f_Rt2H(Rh_m, th_m);
          f_3Dframe(Hh,'b',0.1,'_{h}');
      
%% Image projections on the HELMET cameras
        [u,v]=f_helm_proj(P_p, Hh_m, Hc_h, K, 2);
    pause
%% Compute 3 rays -each one coming from different camera (not concurrent)-       
       % Select three pixels 
        u1 = u(1,1); v1 = v(1,1); %(1st point proj IN 1st camera)
        u2 = u(2,2); v2 = v(2,2); %(2nd point proj IN 2nd camera)
        u3 = u(3,3); v3 = v(3,3); %(3rd point proj IN 3rd camera)
        
       % Unit directions for each pixel (in each camera frame)
        d1_c = inv(K(:,:,1))*[u1;v1;1];
        d1_c = d1_c/norm(d1_c);
        d2_c = inv(K(:,:,2))*[u2;v2;1];
        d2_c = d2_c/norm(d2_c);
        d3_c = inv(K(:,:,3))*[u3;v3;1];
        d3_c = d3_c/norm(d3_c);
       % Unit direction vectors (in <h>)
        d1 = Rc_h(:,:,1)*d1_c;
        d2 = Rc_h(:,:,2)*d2_c;
        d3 = Rc_h(:,:,3)*d3_c;
       % Points on the rays (= tc_h ... for each camera)
        p1 = tc_h(:,1);
        p2 = tc_h(:,2);
        p3 = tc_h(:,3);
        
        
%% 7.1 - Lining-up the rays
    d4 = (f_skew(d1)*d2)/norm(f_skew(d1)*d2);
    d5 = (f_skew(d1)*d4)/norm(f_skew(d1)*d4);
    
    R1 = [ d5 d1 d4 ];
    
    s = d1.'*d2/(d5.'*d2);
        
    alpha = (d1 - s*d5).'*(p2-p1);
    p4 = p1 + alpha*d1;
    
    H1 = [ R1 , -R1*p4;
            [0 0 0], 1];
        
%% 7.2 - Lining-up the Points
     q1 = P_p(:,1);
     q2 = P_p(:,2);
     q3 = P_p(:,3);
     
     e = p2(3);     
     D = sqrt( norm(q2 - q1)^2 - e^2 ); %% IMMAGINARY ??
     
     d6 = [ D; 0; e ]/norm(q2-q1);     
     d7 = (q2 - q1)/norm(q2-q1);     
     d8 = [ 0; 1; 0];     
     d9 = f_skew(d7)*(q3-q1)/norm( f_skew(d7)*(q3-q1) );
     
     R2 = [d6 d8 f_skew(d6)*d8]*[d7 d9 f_skew(d7)*d9];
     
     H2 = [R2  , -R2*q1;
           [ 0 0 0 ], 1];
       
%% 7.3 - Computing the plane coefficients       
    norms= [norm(f_skew(d3)*[1;0;0]), norm(f_skew(d3)*[0;1;0])];
    [max_n, ind_n]= max(norms);
    if ind_n==1,
        n = f_skew(d3)*[1;0;0];
    else
        n = f_skew(d3)*[0;1;0];
    end
    np=f_skew(d3)*n;
    
    L  = [n.'  -n.'*p3]';
    Lp = [np.' -np.'*p3]';
    
%% 7.4 - Intersecting the Quartic and the Circle
     tic
    z=sym('z');
    
     k1 = - e/D*z + norm(q2-q1)*q3'*d6/D;
     k2 = norm(q3)^2 - z^2 - k1^2;
     k3 = L(2);
     k4 = L(1)*k1 + s*D*L(2);
     k5 = -L(1);
     k6 = L(2)*k1 - D*L(2);
     k7 = z*L(3)+L(4);
     k8 = Lp(2);   
     k9 = Lp(1)*k1 + s*D*Lp(2);
    k10 = -Lp(1);
    k11 = Lp(2)*k1 -D*Lp(2);
    k12 = z*Lp(3) + Lp(4);
    k13 = k5*k9 + k6*k8 - k3*k11 - k4*10;
    k14 = k6*k9 - k4*k11 + k2*(k5*k8 - k3*k10);
    k15 = k7*k10 - k5*k12;
    k16 = k7*k11 - k6*k12;
    k17 = k7*k8 - k3*k12;
    k18 = k7*k9 - k4*k12;
    k19 = k2*(k13^2-k15^2-k17^2) + k14^2 - k16^2 - k18^2;
    k20 = 2*(k15*k16 + k17*k18 - k13*k14);
    octic_pol= k19^2  - k2*k20^2;
    coeff_pol=sym2poly(octic_pol);
    beta=coeff_pol/coeff_pol(1);
    
%% 7.5 - Solving the octic
        M = [-beta(2:9);
             1 zeros(1,7);
             0 1 zeros(1,6);
             0 0 1 zeros(1,5);
             0 0 0 1 zeros(1,4);
             0 0 0 0 1 zeros(1,3);
             0 0 0 0 0 1 zeros(1,2);
             0 0 0 0 0 0 1 zeros(1,1)];  % 8 x8 companion matrix
    Z = eig(M);
    time=toc
    
    % Uncomment this to get only solutions related to real coefficients!
    %Zr=Z([5:8]);
    %clear Z
    %Z=Zr;
    
%% 7.6 - Backsubstitution
    for i=1:length(Z);
        X(i) = - e/D*Z(i) + norm(q2-q1)*q3.'*d6/D;
                % Z(i) is varying depending on "i"
                 k1 = - e/D*Z(i) + norm(q2-q1)*q3.'*d6/D;
                 k2 = norm(q3)^2 - Z(i)^2 - k1^2;
                 k3 = L(2);
                 k4 = L(1)*k1 + s*D*L(2);
                 k5 = -L(1);
                 k6 = L(2)*k1 - D*L(2);
                 k7 = Z(i)*L(3)+L(4);
                 k8 = Lp(2);   
                 k9 = Lp(1)*k1 + s*D*Lp(2);
                k10 = -Lp(1);
                k11 = Lp(2)*k1 -D*Lp(2);
                k12 = Z(i)*Lp(3) + Lp(4);
                k13 = k5*k9 + k6*k8 - k3*k11 - k4*10;
                k14 = k6*k9 - k4*k11 + k2*(k5*k8 - k3*k10);
                k15 = k7*k10 - k5*k12;
                k16 = k7*k11 - k6*k12;
                k17 = k7*k8 - k3*k12;
                k18 = k7*k9 - k4*k12;
                k19 = k2*(k13^2-k15^2-k17^2) + k14^2 - k16^2 - k18^2;
                k20 = 2*(k15*k16 + k17*k18 - k13*k14);
                % end
        Y(i) = k19/k20;
        
        
                a(1) = X(i)*L(1) + Y(i)*L(2) + s*D*L(2);            
                a(2) = X(i)*L(2) - Y(i)*L(1) - D*L(2);
                a(3) = Z(i)*L(3) + L(4);
                ap(1) = X(i)*Lp(1) + Y(i)*Lp(2) + s*D*Lp(2);            
                ap(2) = X(i)*Lp(2) - Y(i)*Lp(1) - D*Lp(2);
                ap(3) = Z(i)*Lp(3) + Lp(4);
              u(i) =  (a(2)*ap(1) - ap(2)*a(1))/(ap(2)*a(3)-a(2)*ap(3));
              k(i) = -(a(3)*ap(1) - ap(3)*a(1))/(ap(2)*a(3)-a(2)*ap(3));
            H3([1:4],[1:4],i)=[ 1, -k(i), 0,       0;
                               k(i),  1,   0,   (s-k(i))*D;  
                                 0,   0,  u(i),     0;
                                 0,   0,   0,     u(i)]; %Doubt!!
                         % to check
                          R3([1:3],[1:3],i)= (1/u(i))*[ 1, -k(i), 0;
                                                k(i), 1,   0;  
                                                0,   0,  u(i)];
            d10 = f_skew(d6)*real([X(i);Y(i);Z(i)])/norm(f_skew(d6)*[X(i);Y(i);Z(i)]);  
        
            R4(:,:,i) = [d6 d10 f_skew(d6)*d10]*[d6 d8 f_skew(d6)*d8];
            H4(:,:,i) = [ R4(:,:,i) [0;0;0];
                         [0,0,0], 1 ];
    
         H(:,:,i) = inv(H1)*H3(:,:,i)*H4(:,:,i)*H2;   
         
         
         % Plot the 8 helmet frames
         hold on
          %Hh = f_Rt2H(H([1:3],[1:3],i), H([1:3],4,i));
          f_3Dframe(H(:,:,i),'r',0.2, num2str(i));
  end;
    
           
          
          
    figure,    ezsurf(octic_pol,[-1 1 -1 1])
    
    
    
    
    