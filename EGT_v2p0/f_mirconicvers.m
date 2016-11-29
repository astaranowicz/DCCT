%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  Epipolar Geometry Toolbox  (EGT)  %%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  f_mirconicvers_NEW    Plots the epipolar conic on the mirror.
%  f_mirconicvers_NEW(H,vers,a,b,fig,color)
%
%
%%  Descr: 
%  -----  This function plots, in the 3D space, the epipolar conic 
%         generated on the mirror surface of the current camera, by the 
%         epipolar plane with normal versor VERS.
%
%  H = an homogeneous transformation with respect to the world reference 
%      frame containing both rotation R and translation t, of current 
%      camera. H is a 4by4 matrix (for more details, see the EGT manual).
%  vers = a column vector which defines the normal versor to the epipolar
%         plane, expressed in the mirror frame of the current camera.
%         It is a 3by1 vector (see the EGT manual for more details).
%  a,b = hyperbolic mirror parameters (scalar values).
%  fig = integer to define the figure into which drawing.
%  color = character string made from the characters listed under the 
%          PLOT command.
%
%
%  Example
%     figure(1), hold on
%     axis equal,
%     f_3Dwf(1);
%     H = f_Rt2H(eye(3),[10 0 0]');
%     f_3Dframe(H,1,'g','_d');
%     f_mirconicvers(H,[1 1 1]',3,1,1,'c');
%     f_3Dpanoramic(H,'g',3,1,1,2);


function f_mirconicvers_NEW(H,vers,a,b,fig,color)

if nargin==2
    a = 3;
    b = 1;
    fig = figure;
    color = 'b';
elseif nargin==3
    b = 1;
    fig = figure;
    color = 'b';
elseif nargin==4
    fig = figure;
    color = 'b';
elseif nargin==5
    color = 'b';
elseif nargin>7
    display('  EGT error: too many inputs in f_mirconicvers');
end
    
  
plotta=1;  %togliere 

Rwf2mir = H(1:3,1:3);
twf2mir = H(1:3,4);
vers_wf = Rwf2mir*vers+twf2mir;

r = vers(1);
s = vers(2);
t = vers(3);

% controllo per verificare che 'vers' sia un versore valido
if ((r==0)|(abs(r)<2^-52))&((s==0)|(abs(s)<2^-52))&((t==0)|(abs(t)<2^-52))
    disp(' Il versore seguente non è accettabile come versore normale ad un piano');
    disp(vers);
else
    if (r==0)|(abs(r)<2^-52)     % se r è uguale a 0
        if (s>=0)|(abs(s)<2^-52)
            alfa = 0;
        else  % se t<0
            alfa = pi;
        end
    else                         % se r è diverso da 0
        if (s==0)|(abs(s)<2^-52)
            if r>0
                alfa = pi/2;
            else  % if r<0
                alfa = -pi/2;
            end   
        elseif s>0
            alfa = atan(r/s);
        else  % se t<0
            alfa = pi+atan(r/s);
        end
    end   
    
    if (t==0)|(abs(t)<2^-52)
        beta = pi/2;
    elseif t>0
        beta = atan(sqrt(r^2+s^2)/t);
    else  % se s<0
        beta = pi+atan(sqrt(r^2+s^2)/t);
    end
end
        
A = rotoy(alfa);
B = rotox(beta);
R_ab = A*B;
P = [R_ab      , zeros(3,1); 
     zeros(1,3)      1     ];
e = sqrt(a^2+b^2);

% % if t==0
% %     if r>0
% %         phi = pi/2;
% %     elseif r<0
% %         phi = (3/2)*pi;
% %     else
% %         phi = 0;
% %     end
% % else
% %     phi = atan(r/t);
% % end;
% % if s==0
% %     if t>=0
% %         psi = pi/2;
% %     else
% %         psi = -pi/2;
% %     end
% % else
% %     psi = atan(sqrt(r^2+t^2)/s);
% % end;
% % A = rotoy(phi);
% % B = rotox(psi);
% % R01 = A*B;
% % P = [R01, zeros(3,1); zeros(1,3) 1];
% % e = sqrt(a^2+b^2);

%

%%%NEW
Q0=[ -a^2    0       0     0
       0    -a^2     0     0
       0     0      b^2   e*b^2
       0      0    e*b^2   b^4];       
                            
% Q0 = [ -a^2    0    0    0   ;
%          0    b^2   0  e*b^2 ;
%          0     0  -a^2   0   ;
%          0   e*b^2  0   b^4  ];

Q1 = P'*Q0*P;
%Q1
%pause
% C2D = [Q1(3,3) Q1(3,1) Q1(3,4);
%        Q1(1,3) Q1(1,1) Q1(1,4);
%        Q1(4,3) Q1(4,1) Q1(4,4)];
 
C2D = [Q1(1,1) Q1(1,3) Q1(1,4);
       Q1(3,1) Q1(3,3) Q1(3,4);
       Q1(4,1) Q1(4,3) Q1(4,4)];

detC2D = det(C2D);
delta = det(C2D(1:2,1:2));
den = (C2D(2,2)-C2D(1,1));
if den==0
    phi = pi/2;
else
    phi = (atan(2*C2D(1,2)/(C2D(2,2)-C2D(1,1))))/2;
end;
% u = R01*u1  ==> u1 = R01'*u 
rotphi = rotoz(phi);%era "z"
R01 = rotphi;
C2D1 = rotphi'*C2D*rotphi;

if detC2D~=0
    if delta~=0
        u01 = -(C2D1(1,3)/C2D1(1,1));
        v01 = -(C2D1(2,3)/C2D1(2,2));
    else
        if C2D1(2,2)==0
            u01 = -(C2D1(1,3)/C2D1(1,1));
            v01 = (-C2D1(1,1)*u01^2-2*C2D1(1,3)*u01-C2D1(3,3))/(2*C2D1(2,3));
        elseif C2D1(1,1)==0
            v01 = -(C2D1(2,3)/C2D1(2,2));
            u01 = (-C2D1(2,2)*v01^2-2*C2D1(2,3)*v01-C2D1(3,3))/(2*C2D1(1,3));
        end
    end
    
    C2D2 = [        C2D1(1,1)                    0                    C2D1(1,1)*u01+C2D1(1,3);
                        0                    C2D1(2,2)                C2D1(2,2)*v01+C2D1(2,3);
            C2D1(1,1)*u01+C2D1(1,3)  C2D1(2,2)*v01+C2D1(2,3)  C2D1(1,1)*u01^2+C2D1(2,2)*v01^2+C2D1(3,3)+2*C2D1(1,3)*u01+2*C2D1(2,3)*v01];
    
    if delta>0
        [u2,v2] = f_ellipse(C2D2);
        u1 = u2+u01;
        v1 = v2+v01;
        u = R01(1,:)*[u1;v1;ones(1,length(u1))];
        v = R01(2,:)*[u1;v1;ones(1,length(u1))]; 
        % [x,0,z,1]' ---> [u,0,v,1]'
        ellisse = [Rwf2mir*(A*B) twf2mir; zeros(1,3) 1]*[u;zeros(1,length(v));v;ones(1,length(v))]; 
        if plotta==1
          figure(fig),plot3(ellisse(1,:),ellisse(2,:),ellisse(3,:),color); 
        end
    elseif delta<0 
        [u2,v2] = f_hyperbola(C2D2);
        u1 = u2+u01;
        v1 = v2+v01;
        u = R01(1,:)*[u1;v1;ones(1,length(u1))];
        v = R01(2,:)*[u1;v1;ones(1,length(u1))];
        u1neg = -u2+u01; 
        uneg = R01(1,:)*[u1neg;v1;ones(1,length(u1))]; 
        vneg = R01(2,:)*[u1neg;v1;ones(1,length(u1))];
        iperbole = [Rwf2mir*(A*B) twf2mir; zeros(1,3) 1]*[u;zeros(1,length(v));v;ones(1,length(v))];  
        iperbole1 = [Rwf2mir*(A*B) twf2mir; zeros(1,3) 1]*[uneg;zeros(1,length(v));vneg;ones(1,length(v))];  
        if plotta==1
            figure(fig),hold on,plot3(iperbole(1,:),iperbole(2,:),iperbole(3,:),color);
            plot3(iperbole1(1,:),iperbole1(2,:),iperbole1(3,:),color);
        end
    elseif delta==0
        [u2,v2] = f_parabola(C2D2);
        u1 = u2+u01;
        v1 = v2+v01;
        u = R01(1,:)*[u1;v1;ones(1,length(u1))];
        v = R01(2,:)*[u1;v1;ones(1,length(u1))];
        parabola = [Rwf2mir*(A*B) twf2mir; zeros(1,3) 1]*[u;zeros(1,length(v));v;ones(1,length(v))]; 
        if plotta==1
            figure(fig),plot3(parabola(1,:),parabola(2,:),parabola(3,:),color);       
        end
    end
%     figure(fig), hold on
%     plot3([twf2mir(1) vers_wf(1)],[twf2mir(2) vers_wf(2)],[twf2mir(3) vers_wf(3)],'k')
%     col = strcat('k','>');
%     plot3(vers_wf(1), vers_wf(2), vers_wf(3),col)
else
    display('  conica degenere');  % coniche degeneri
end;   
