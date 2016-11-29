function f_3Dcamera_helm(Hh_m, Hc_h, scale, color);
%Camera points in the <wf>
P_c=[0,  0, 0;
     1, -1, 2;
     1,  1, 2;
    -1,  1, 2;
    -1, -1, 2]*scale;
n=length(P_c(:,1));

% These are the cameras 3D points belonging to  the head, expressed in the
% head frame...
    Rc_h = Hc_h([1:3],[1:3]);
    tc_h = Hc_h([1:3],4);
    Tc_h = tc_h*ones(1,n);
P_h= (Rc_h*P_c' + Tc_h)';

% for the plot I have to express them back in the world frame
        Rh_w = Hh_m([1:3],[1:3]); %this is given as a rotation from <w> to <h> in reality
        Rw_m=rotox(-pi/2);
    Rh_m = Rw_m*Rh_w;
    th_m = Hh_m([1:3],4); %this instead in centered in the <M> frame
    Th_m = th_m*ones(1,n);
P=(Rh_m*P_h' + Th_m)';

%COmponents of the 1st patch
X1=[P(2,1) ; P(3,1) ; P(4,1) ; P(5,1)];
Y1=[P(2,2) ; P(3,2) ; P(4,2) ; P(5,2)];
Z1=[P(2,3) ; P(3,3) ; P(4,3) ; P(5,3)];
%COmponents of the 2nd patch
X2=[P(1,1) ; P(2,1) ; P(3,1) ; P(1,1)];
Y2=[P(1,2) ; P(2,2) ; P(3,2) ; P(1,2)];
Z2=[P(1,3) ; P(2,3) ; P(3,3) ; P(1,3)];
%COmponents of the 2nd patch
X3=[P(1,1) ; P(3,1) ; P(4,1) ; P(1,1)];
Y3=[P(1,2) ; P(3,2) ; P(4,2) ; P(1,2)];
Z3=[P(1,3) ; P(3,3) ; P(4,3) ; P(1,3)];
%COmponents of the 2nd patch
X4=[P(1,1) ; P(4,1) ; P(5,1) ; P(1,1)];
Y4=[P(1,2) ; P(4,2) ; P(5,2) ; P(1,2)];
Z4=[P(1,3) ; P(4,3) ; P(5,3) ; P(1,3)];
%COmponents of the 2nd patch
X5=[P(1,1) ; P(5,1) ; P(2,1) ; P(1,1)];
Y5=[P(1,2) ; P(5,2) ; P(2,2) ; P(1,2)];
Z5=[P(1,3) ; P(5,3) ; P(2,3) ; P(1,3)];

hold on;
color=[.2 .2 .2];
opaqueness=.3;
fill3(X1,Y1,Z1,[.9 .9 .9],'FaceAlpha',.8);
fill3(X2,Y2,Z2,color+[.5 .5 .5],'FaceAlpha',opaqueness);
fill3(X3,Y3,Z3,color+[.3 .3 .3],'FaceAlpha',opaqueness);
fill3(X4,Y4,Z4,color+[.1 .1 .1],'FaceAlpha',opaqueness);
fill3(X5,Y5,Z5,color,'FaceAlpha',opaqueness);
hold off;