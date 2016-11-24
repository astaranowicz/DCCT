function [u,v]=f_helm_proj(P_p, Hh_m, Hc_h, K, b_plot);

        Rh_m=Hh_m([1:3],[1:3]);
        th_m=Hh_m([1:3],4);
        Rh_w=Rh_m;
        Rm_w=rotox(pi/2);%ok
        th_w=Rm_w*th_m;%ok
        
        
        for i=1:length(Hc_h(1,1,:)),
            tc_h(:,i)=Hc_h([1:3],4,i);
            Rc_h(:,:,i)=Hc_h([1:3],[1:3],i);
            tc_w(:,i) = Rh_w*tc_h(:,i)+th_w;%ok
            Rc_w(:,:,i) = Rh_w*Rc_h(:,:,i);
            tc_m(:,i) = Rm_w'*tc_w(:,i);
            Hc([1:4],[1:4],i)= [  Rc_w(:,:,i) , tc_m(:,i);
                                   [ 0 0 0] ,   1     ];  
            [u(i,:),v(i,:)]=f_perspproj( P_p, Hc(:,:,i), K(:,:,i), b_plot);
        end;