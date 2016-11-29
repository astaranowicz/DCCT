%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C) 2016 Aaron Staranowicz and Gian Luca Mariottini
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% LEAST-SQUARES LINE FITTING
%
% [a,b]=f_linefitLS(x,y,alg);
% where
%   [a b] = line parameters (ax+b=0)
%   (x,y) = 2-D point set (with noise)
%   alg=1 : Pseudo inverse (minimization of orthogonal distances)
%   alg=2 : minimization of vertical distances
%
% INFO: See my notes on line fitting at pag. B
%
function [a,b]=f_linefitLS(x,y,alg);
m=length(x);
if alg==1,
    A(:,[1:2]) = [ones(m,1) x];
  	Y(:,1)=y;
	%X_hat=inv(A'*A)*A'*Y;
    X_hat = A\Y;
    b=X_hat(1);
    a=X_hat(2);
elseif alg==2,
    Av = [ sum(x.^2) sum(x);
            sum(x)   length(x) ];
	  Yv = [ sum(y.*x); sum(y) ];
    %Xv_hat=inv(Av)*Yv;
    Xv_hat = Av\Yv;
    a=Xv_hat(1);
    b=Xv_hat(2);
elseif alg==3,
    %% ALG3: Line fitting with m points (m>=3) for the
    %% (cos(theta),sin(theta),rho) parameterization
    if m<3,
        display('ERROR "f_linefitLS" : Alg.3 needs at least three points')
        return
    else
        A(:,[1:3]) = [x y ones(m,1)];
        [Uunr,Dunr,Vunr]=svd(A'*A);
        Dunr(3,3)=0; %force the rank 2 constraint
        At=Uunr*Dunr*Vunr';
        [U,D,V]=svd(At);
        theta_hat=V(:,end);
        theta_hat=theta_hat*sign(theta_hat(3));
        scale=sqrt(theta_hat(1)^2+theta_hat(2)^2); %Bring the third component to be always positive
        theta_hat=theta_hat/scale;

    % Obtain (a,b) of the line from (rho,theta) 
        cos_t=theta_hat(1);
        sin_t=theta_hat(2);
        rho  =theta_hat(3);
        theta= atan2(sin_t,cos_t);
        theta_1=asin(sin_t);
        theta_2=acos(cos_t);
        theta=(sign(sin_t)*sign(cos_t))*(abs(theta_1)+abs(theta_2))/2;
        
        % Convert from [rho,theta] to [a,b]
        v=[rho*cos(theta); rho*sin(theta)];
        d=[-sin(theta); +cos(theta)];
        P1=v;
        P2=v+d;
        % Find (a,b) for the line passign through P1 and P2
        A =[P1(1)^2+P2(1)^2,   P1(1)+P2(1);
              P1(1)+P2(1),          2    ];
        z =[P1(1)*P1(2)+P2(1)*P2(2);
            P1(2)+P2(2)];
        a_b = A\z;%inv(A)*z;
        a=a_b(1);
        b=a_b(2);
    %% However at the end there might be an ambiguity on theta of +/- pi, easily solved
    %% counting the inliers
        T=1; %Threshold for checking how many inliers
        dist=abs(a*x+b-y); %vector of distances
        if length(find(dist<=T)) == 0
            % Then invert the sign, which corresponds on taking the other
            % angle (theta+pi)
            v=[rho*(-cos(theta)); rho*(-sin(theta))];
            d=[-(-sin(theta)); (-cos(theta))];
            P1=v;
            P2=v+d;
            % Find (a,b) for the line passign through P1 and P2
            A =[P1(1)^2+P2(1)^2,   P1(1)+P2(1);
                  P1(1)+P2(1),          2    ];
            z =[P1(1)*P1(2)+P2(1)*P2(2);
                P1(2)+P2(2)];
            a_b = A\z;%inv(A)*z;
            a=a_b(1);
            b=a_b(2);
        end
    end
%% Minimal polar line parameteriz. fitting
elseif alg==4,    
    %rho_meas=sqrt(x.^2+y.^2);
    %theta_meas=atan2(y,x);
    
    [theta_meas,rho_meas]=cart2pol(x,y);
    
    A=zeros(m,2);
    for i=1:m,
        A(i,[1:2])=rho_meas(i)*[cos(theta_meas(i)), sin(theta_meas(i))];
    end
    % A.x=1
    sigma=1;
    W=1/sigma^2*diag(m);
    %x=inv(A'*A)*A'*ones(m,1);
    Xx=A\ones(m,1)
    
    rho=sqrt(Xx(1)^2+Xx(2)^2)
    theta=atan2(Xx(2),Xx(1))
      
    % Convert from [rho,theta] to [a,b]
        v=[rho*cos(theta); rho*sin(theta)];
        d=[-sin(theta); +cos(theta)];
        P1=v;
        P2=v+d;
        % Find (a,b) for the line passign through P1 and P2
        Aa =[P1(1)^2+P2(1)^2,   P1(1)+P2(1);
              P1(1)+P2(1),          2    ];
        z =[P1(1)*P1(2)+P2(1)*P2(2);
            P1(2)+P2(2)];
        a_b = Aa\z;%inv(A)*z;
        a=a_b(1);
        b=a_b(2);
    
end
