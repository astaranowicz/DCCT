%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Epipolar Geometry Toolbox v1.3 (EGT)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [Xcd_e,e_cd]=f_panepipoles(Hd,Ha,a1,b1,quadric1,a2,b2,quadric2,K1,K2,color1,color2);
%
%- DESCRIPTION:
%         "f_panepipoles" computes the mirror (<wf>) and image (pixels) coordinates of the epipoles 
%         existing between two central catadioptric cameras. Mirror epipoles are also plotted.
%
%- SYNTAX:
%  Xcd_e = first output 3x4 matrix containing the four epipoles in the <w.f.>
%         (onto the mirror surface) (note that Xw_d1 is a 3x1 vector):
%                   Xcd_e = [ Xw_d1 , Xw_d2 , Xw_c1 , Xw_c2 ];
%  e_cd  = second output 3x4 matrix containing the four epipoles in the
%          image plane (note that e_d1 is a 2x1 vector):
%                    e_cd = [ e_d1 , e_d2 , e_c1 , e_c2 ];
%  Hd, Ha = Homogeneous matrix containing rotation R and translation t,
%           defined with respect to the world frame, of desired and 
%           actual camera respectively (for more details see the EGT manual). 
%  a1, b1, quadric1 (a2, b2, quadric2) 
%         = mirror parameters of desired and current camera respectively. 
%  K1, K2 = calibration matrix with the inner parameters of desired
%           and actual camera. K1 and K2 are 3x3 matrices.
%  color1, color2 = character string which defines the colour used to draw
%                   epipoles of desired and actual camera respectively
%                   (e.g., 'g', 'b', etc.)
%   
%- USAGE:
%  [e1d,e2d,e1a,e2a] = f_panepipoles(Hd,Ha) returns the pixel coordinates
%  of the epipoles of two panoramic cameras positioned how defined in Hd 
%  e Ha respectively and with the other parameters having the default value: 
%              a1 = a2 = 3, b1 = b2 = 1, Kd = Ka = eye(3).
%
%Authors:
%  Gian Luca Mariottini, Eleonora Alunno
%
%Last update:
%  December 2005
%

function [Xcd_e,e_cd]=f_panepipoles(Hd,Ha,a1,b1,quadric1,a2,b2,quadric2,K1,K2,color1,color2);
if nargin<2
    display('EGT error: f_panepipoles has insufficient inputs!')
elseif nargin==2,
    a1=0.04; b1=0.02; quadric1=1; a2=0.04; b2=0.02; quadric2=1;
    K1=[10^3  0   640; 0  10^3  480;  0    0    1 ]; K2=[10^3  0   640; 0  10^3  480;  0    0    1 ];
    color1='g'; color2='r';
elseif nargin==3
    b1=0.02; quadric1=1; a2=0.04; b2=0.02; quadric2=1;
    K1=[10^3  0   640; 0  10^3  480;  0    0    1 ]; K2=[10^3  0   640; 0  10^3  480;  0    0    1 ];
    color1='g'; color2='r';
elseif nargin==4,
    quadric1=1; a2=0.04; b2=0.02; quadric2=1;
    K1=[10^3  0   640; 0  10^3  480;  0    0    1 ]; K2=[10^3  0   640; 0  10^3  480;  0    0    1 ];
    color1='g'; color2='r';
elseif nargin==5,
    a2=0.04; b2=0.02; quadric2=1;
    K1=[10^3  0   640; 0  10^3  480;  0    0    1 ]; K2=[10^3  0   640; 0  10^3  480;  0    0    1 ];
    color1='g'; color2='r';
elseif nargin==6,
    b2=0.02; quadric2=1;
    K1=[10^3  0   640; 0  10^3  480;  0    0    1 ]; K2=[10^3  0   640; 0  10^3  480;  0    0    1 ];
    color1='g'; color2='r';
elseif nargin==7,
    quadric2=1;
    K1=[10^3  0   640; 0  10^3  480;  0    0    1 ]; K2=[10^3  0   640; 0  10^3  480;  0    0    1 ];
    color1='g'; color2='r'; 
elseif nargin==8,
    K1=[10^3  0   640; 0  10^3  480;  0    0    1 ]; K2=[10^3  0   640; 0  10^3  480;  0    0    1 ];
    color1='g'; color2='r'; 
elseif nargin==9,
    K2=[10^3  0   640; 0  10^3  480;  0    0    1 ];
    color1='g'; color2='r'; 
elseif nargin==10,
    color1='g'; color2='r'; 
elseif nargin==11,
    color2='r'; 
elseif nargin>12,
    display('EGT error: f_panepipoles has more outputs than expected!');
end;    
    % world frame == matlab frame
    td_w  = Hd([1:3],4); % td_w stands for "from <w> to <d>"
    tc_w  = Ha([1:3],4);      
    Rd_w = Hd([1:3],[1:3]);
    Rc_w = Ha([1:3],[1:3]);   
    Rm_c = eye(3);
    tm_c1 = [0 0 2*sqrt(a1^2+b1^2)]';
    tm_c2 = [0 0 2*sqrt(a2^2+b2^2)]';
    Rw_d = Rd_w';
    Rw_c = Rc_w';
    tdc_w= tc_w - td_w; %From <d> to <c> as it is in the world frame
    tc_d = Rw_d*tdc_w;
    td_c = Rw_c*(-tdc_w);

    [lambda_d1,lambda_d2] = f_lambda(a1,b1,quadric1,tc_d);
    [lambda_c1,lambda_c2] = f_lambda(a2,b2,quadric2,td_c);
    
    %% Epipoles in the mirror surface expressed in the <d> and <c> frames
    X_d1=lambda_d1*tc_d;
    X_d2=lambda_d2*tc_d;
    X_c1=lambda_c1*td_c;
    X_c2=lambda_c2*td_c;
    
    %% Projection of epipoles to the image plane (orthographic or
    %% perspective, wrt the value of quadric).
    % Formulas are taken from Mariottini, 2005. Ph.D. Thesis
    if quadric1==2,        
        e_d1=K1*[ [1 0 0; 0 1 0]*X_d1;
                       1      ]; %Image plane projection of epipole originated by X_d1 (mirr. frame)
        e_c1=K1*[ [1 0 0; 0 1 0]*X_c1;
                       1      ];                              
    end
    if quadric2==2,
        e_d2=K2*[ [1 0 0; 0 1 0]*X_d2;
                       1      ];
        e_c2=K2*[ [1 0 0; 0 1 0]*X_c2;
                       1      ];           
    end;
    if quadric1==1,
        Xc1=(Rm_c*X_d1+tm_c1);
        e_d1=K1*(1/Xc1(3))*(X_d1+[0;0;2*sqrt(a1^2+b1^2)]);
        e_c1=K1*(1/Xc1(3))*(X_c1+[0;0;2*sqrt(a1^2+b1^2)]);
    end    
    if quadric2==1,
        Xc2=(Rm_c*X_d2+tm_c2);
        e_d2=K2*(1/Xc2(3))*(X_d2+[0;0;2*sqrt(a2^2+b2^2)]);
        e_c2=K2*(1/Xc2(3))*(X_c2+[0;0;2*sqrt(a2^2+b2^2)]);
    end    
    
    %% Express the epipoles in the world frame (onto the mirror surface)
    Xw_c1=Rc_w*X_c1+tc_w;
    Xw_c2=Rc_w*X_c2+tc_w;
    Xw_d1=Rd_w*X_d1+td_w;
    Xw_d2=Rd_w*X_d2+td_w;
    
    %% OUTPUTS    
    Xcd_e=[Xw_d1,Xw_d2,Xw_c1,Xw_c2];
    e_cd =[e_d1 ,e_d2 ,e_c1 ,e_c2 ];

    %% Plot epipoles in the world frame
    plot3(Xw_c1(1),Xw_c1(2),Xw_c1(3),strcat(color2,'*'));
    plot3(Xw_c2(1),Xw_c2(2),Xw_c2(3),strcat(color2,'*'));
    plot3(Xw_d1(1),Xw_d1(2),Xw_d1(3),strcat(color1,'*'));
    plot3(Xw_d2(1),Xw_d2(2),Xw_d2(3),strcat(color1,'*'));   
    text(Xw_c1(1),Xw_c1(2),Xw_c1(3),'e_{c1}');
    text(Xw_c2(1),Xw_c2(2),Xw_c2(3),'e_{c2}');
    text(Xw_d1(1),Xw_d1(2),Xw_d1(3),'e_{d1}');
    text(Xw_d2(1),Xw_d2(2),Xw_d2(3),'e_{d2}');
    
    plot3([Xw_c1(1) Xw_d1(1)],[Xw_c1(2) Xw_d1(2)],[Xw_c1(3) Xw_d1(3)],'m');