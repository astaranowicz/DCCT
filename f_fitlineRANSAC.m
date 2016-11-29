function [m,q,inliers] = f_fitlineRANSAC(x,y,thresh,alg)

n=length(x);
npunti=2;
list_points = nchoosek([1:n],npunti);
N=nchoosek(n,npunti);
consensus=zeros(1,N);
for i=1:N,
    p1_x = x(list_points(i,1));   
    p1_y = y(list_points(i,1));    
    p2_x = x(list_points(i,2)); 
    p2_y = y(list_points(i,2));    
    %keyboard
    [m1,q1] = f_linefitLS([p1_x p2_x]', [p1_y p2_y]', alg);    
    
    for j=1:n,
        theta=pi/2+atan(m1);
        rho=q1*sin(theta);
        %dist(j) = (y(j)-(m2*x(j)+q2))^2;
        dist(j)= abs(x(j)*cos(theta)+y(j)*sin(theta)-rho);
        if dist(j)<= thresh,
          consensus(i) = consensus(i)+1;
        end
    end 

    %dist
    %pause
    %consensus(i)
    %pause
end


% Take the one with highest consensus
[max_cons, ind_max_1] = max(consensus);
% Find if there are other consensus with the same cost
inds_max = find(consensus==max_cons);

    for p=1:length(inds_max),
        ind_max = inds_max(1);
        % Use the max to estimate again and take inliers
        p1_x = x(list_points(ind_max,1));
        p1_y = y(list_points(ind_max,1));
        p2_x = x(list_points(ind_max,2));    
        p2_y = y(list_points(ind_max,2));
        [m2,q2] = f_linefitLS([p1_x p2_x]', [p1_y p2_y]',alg);
        inliers=[];
        for j=1:n,
            theta=pi/2+atan(m2);
            rho=q2*sin(theta);
            %dist(j) = (y(j)-(m2*x(j)+q2))^2;
            dist(j)= abs(x(j)*cos(theta)+y(j)*sin(theta)-rho);
            if dist(j)<= thresh,  
              inliers = [inliers, j ];
            end
        end 
        % Final estimate
        [m,q] = f_linefitLS(x(inliers), y(inliers), alg);
%         plot([2 8], [2*m+q, 8*m+q],'g');
    end
    
%     figure(100)
%       clf(100)
%     plot([0 6], [q, 6*m+q],'r');
%     hold on
%     axis equal
%     plot(x,y,'b.')
%     plot(x(inliers),y(inliers),'gO')
%     pause
    
end




  
