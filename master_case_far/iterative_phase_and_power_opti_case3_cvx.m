function [v, p, status_cvx] = iterative_phase_and_power_opti_case3_cvx(w, r_min, PS1, PS2, PB1, PB2,...
            ch_BS1_to_EU, ch_BS2_to_EU, ch_RIS3_to_EU, ch_BS1_to_RIS3, ch_BS2_to_RIS3, ch_SU1_to_EU, ...
            ch_SU2_to_EU, ch_BS1_to_SU1,  ch_BS2_to_SU2, ch_RIS3_to_SU1, ch_SU1_to_RIS3, ch_RIS3_to_SU2, ch_SU2_to_RIS3, noise_power, nIterations)
%         noise_power = 1;

[N,M] = size(ch_BS1_to_RIS3);
Theta_mat_BS1_RIS_SU1 = diag(ch_RIS3_to_SU1')*ch_BS1_to_RIS3;
Theta_mat_BS1_RIS_EU = diag(ch_RIS3_to_EU')*ch_BS1_to_RIS3;
theta_vec_SU1_RIS_EU = diag(ch_RIS3_to_EU')*ch_SU1_to_RIS3;

Theta_mat_BS2_RIS_SU2 = diag(ch_RIS3_to_SU2')*ch_BS2_to_RIS3;
Theta_mat_BS2_RIS_EU = diag(ch_RIS3_to_EU')*ch_BS2_to_RIS3;
theta_vec_SU2_RIS_EU = diag(ch_RIS3_to_EU')*ch_SU2_to_RIS3;

r_min_EU = 2^r_min -1;

t_next = 0.05*ones(2,1);
p_next = 0.00005*ones(2,1);
eta_next = 0.005*ones(3,1);
pt_next = 0.000005*ones(8,1);
qt_next = 0.000005*ones(8,1);
x_next = 0.05*ones(N,1);
y_next = 0.05*ones(N,1);


tol = 1e-2;     
% cvx_solver Mosek
% cvx_solver sedumi
cvx_solver sdpt3
cvx_precision high
% cvx_solver_settings('MSK_IPAR_INFEAS_REPORT_AUTO','MSK_ON','MSK_IPAR_INFEAS_REPORT_LEVEL',1)
while(1)
    sumrate = zeros(nIterations,1);
    seqobj = zeros(nIterations,1);
    for iIteration=1:nIterations % for SCA algorithm
        cvx_begin quiet  
        variables   x(N,1) y(N,1) z t(2,1) eta(3,1) pt(8,1) qt(8,1) p(2,1) u;
        variable    v(N,1) complex;

        maximize( z + sum( x_next.^2 + y_next.^2 + 2.*x_next.*(x-x_next) + 2*y_next.*(y-y_next)));
         
        subject to
        % Objective function
        norm([t(1)-t(2); 2*z]) <= t(1)+t(2);
        
        for i=1:N
            x(i) == real(v(i));
            y(i) == imag(v(i));
        end
        % UE1 Eq. 42f      
        norm([sqrt(noise_power); ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1)]) <= sqrt(eta(1));    
        
        % UE1 Eq. 42g      
        norm([sqrt(noise_power); ch_BS2_to_SU2'*w(:,2) + v'*Theta_mat_BS2_RIS_SU2*w(:,2)]) <= sqrt(eta(2));
        
        % Edge UE Eq. 42j
        norm([sqrt(noise_power); ch_BS2_to_EU'*w(:,2)+ v'*Theta_mat_BS2_RIS_EU*w(:,2); ch_BS1_to_EU'*w(:,1) + v'*Theta_mat_BS1_RIS_EU*w(:,1)]) <= sqrt(eta(3)); 
        
        % Eq. 42h
        pt(7) == real(v'*theta_vec_SU1_RIS_EU);
        qt(7) == imag(v'*theta_vec_SU1_RIS_EU);
        pt(8) == real(v'*theta_vec_SU2_RIS_EU);
        qt(8) == imag(v'*theta_vec_SU2_RIS_EU);
        (p(1)*norm(ch_SU1_to_EU)^2 + p(1)*(pt_next(7)^2 + qt_next(7)^2) + 2*p_next(1)*pt_next(7)*(pt(7) - pt_next(7)) +2*p_next(1)*qt_next(7)*(qt(7) - qt_next(7))...
       + p(2)*norm(ch_SU2_to_EU)^2 + p(2)*(pt_next(8)^2 + qt_next(8)^2) + 2*p_next(2)*pt_next(8)*(pt(8) - pt_next(8)) +2*p_next(2)*qt_next(8)*(qt(8) - qt_next(8)) + noise_power*u) >= (noise_power*r_min_EU);
       
       % Eq. 42k and 42l
        for i=1:N
            norm([x(i) ; y(i)]) <=1;
        end
        % Eq. 42m & 42n 
        p(1) <= PS1; 
        p(2) <= PS2; 
        

        % Eq. 42b
        pt(1) == real(ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1));
        qt(1) == imag(ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1));
        1e15*(pt_next(1)*(2*pt(1)-pt_next(1)) + qt_next(1)*(2*qt(1) - qt_next(1))) >= 1e15*((t(1)-1)*noise_power);
%         2*pt_next(1)*(pt(1)-pt_next(1))/t_next(1) + 2*qt_next(1)*(qt(1) - qt_next(1))/t_next(1) + ((pt_next(1)^2 + qt_next(1)^2)/t_next(1))*(1-(t(1) - t_next(1))/t_next(1))>= noise_power;
        
        % Eq. 42c
        pt(2) == real(ch_BS2_to_SU2'*w(:,2) + v'*Theta_mat_BS2_RIS_SU2*w(:,2));
        qt(2) == imag(ch_BS2_to_SU2'*w(:,2) + v'*Theta_mat_BS2_RIS_SU2*w(:,2));
        1e15*(pt_next(2)*(2*pt(2)-pt_next(2)) + qt_next(2)*(2*qt(2) - qt_next(2))) >= 1e15*((t(2)-1)*noise_power);
%         2*pt_next(2)*(pt(2)-pt_next(2))/t_next(2) + 2*qt_next(2)*(qt(2) - qt_next(2))/t_next(2) + ((pt_next(2)^2 + qt_next(2)^2)/t_next(2))*(1-(t(2) - t_next(2))/t_next(2))>= noise_power;
       
%         % Eq. 42d
        pt(3) == real(ch_BS1_to_SU1'*w(:,3) + v'*Theta_mat_BS1_RIS_SU1*w(:,3));
        qt(3) == imag(ch_BS1_to_SU1'*w(:,3) + v'*Theta_mat_BS1_RIS_SU1*w(:,3));
        pt_next(3)*(2*pt(3)-pt_next(3)) + qt_next(3)*(2*qt(3) - qt_next(3)) >= eta(1)*r_min_EU ;    
%         2*pt_next(3)*(pt(3)-pt_next(3))/eta_next(1) + 2*qt_next(3)*(qt(3) - qt_next(3))/eta_next(1) + ((pt_next(3)^2 + qt_next(3)^2)/eta_next(1))*(1-(eta(1) - eta_next(1))/eta_next(1))>= r_min_EU;

%         % Eq. 42e
        pt(4) == real(ch_BS2_to_SU2'*w(:,4) + v'*Theta_mat_BS2_RIS_SU2*w(:,4));
        qt(4) == imag(ch_BS2_to_SU2'*w(:,4) + v'*Theta_mat_BS2_RIS_SU2*w(:,4));
        pt_next(4)*(2*pt(4)-pt_next(4)) + qt_next(4)*(2*qt(4) - qt_next(4)) >= eta(2)*r_min_EU ;  
%         2*pt_next(4)*(pt(4)-pt_next(4))/eta_next(2) + 2*qt_next(4)*(qt(4) - qt_next(4))/eta_next(2) + ((pt_next(4)^2 + qt_next(4)^2)/eta_next(2))*(1-(eta(2) - eta_next(2))/eta_next(2))>= r_min_EU;        
     
        % Eq. 42i
        pt(5) == real(ch_BS1_to_EU'*w(:,3) + v'*Theta_mat_BS1_RIS_EU*w(:,3));
        qt(5) == imag(ch_BS1_to_EU'*w(:,3) + v'*Theta_mat_BS1_RIS_EU*w(:,3));
        pt(6) == real(ch_BS2_to_EU'*w(:,4) + v'*Theta_mat_BS2_RIS_EU*w(:,4));
        qt(6) == imag(ch_BS2_to_EU'*w(:,4) + v'*Theta_mat_BS2_RIS_EU*w(:,4));
        (2*pt_next(5)*(pt(5)-pt_next(5))/eta_next(3) + 2*qt_next(5)*(qt(5) - qt_next(5))/eta_next(3) + ((pt_next(5)^2 + qt_next(5)^2)/eta_next(3))*(1-(eta(3) - eta_next(3))/eta_next(3))+ 2*pt_next(6)*(pt(6)-pt_next(6))/eta_next(3) + 2*qt_next(6)*(qt(6) - qt_next(6))/eta_next(3) + ((pt_next(6)^2 + qt_next(6)^2)/eta_next(3))*(1-(eta(3) - eta_next(3))/eta_next(3))) >= u ;
        

        p(1) >= 0;
        p(2) >= 0;
        t(1) >=  1;
        t(2) >=  1;
        u >= 0;
        eta(1) >= 0;
        eta(2) >= 0;
        eta(3) >= 0;
        z>=0;
  
        cvx_end;
        if(strfind(cvx_status,'Solved'))
            t_next = double(t);
            p_next = p;
            eta_next(3) = double(eta(3));
            pt_next = double(pt);
            qt_next = double(qt);
            x_next = double(x);
            y_next = double(y);

            
            norm(w(:,1))^2 + norm(w(:,3))^2;
            norm(w(:,2))^2 + norm(w(:,4))^2;
            % compute the Edge user's rate for each iteration
            SU1_rate = log2(1 + norm(ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1))^2/noise_power);
            SU2_rate = log2(1 + norm(ch_BS2_to_SU2'*w(:,2) + v'*Theta_mat_BS2_RIS_SU2*w(:,2))^2/noise_power);
            % compute the Edge user's rate for each iteration
            EU_rate_at_SU1 = log2(1 + norm(ch_BS1_to_SU1'*w(:,3) + v'*Theta_mat_BS1_RIS_SU1*w(:,3))^2/(norm(ch_BS1_to_SU1'*w(:,1)+ v'*Theta_mat_BS1_RIS_SU1*w(:,1))^2 + noise_power));
            EU_rate_at_SU2 = log2(1 + norm(ch_BS2_to_SU2'*w(:,4) + v'*Theta_mat_BS2_RIS_SU2*w(:,4))^2/(norm(ch_BS2_to_SU2'*w(:,2)+ v'*Theta_mat_BS2_RIS_SU2*w(:,2))^2 + noise_power));
            EU_rate = log2(1+ (norm(ch_BS1_to_EU'*w(:,3) + v'*Theta_mat_BS1_RIS_EU*w(:,3))^2 + ...
                norm(ch_BS2_to_EU'*w(:,4)+ v'*Theta_mat_BS2_RIS_EU*w(:,4))^2)/(norm(ch_BS1_to_EU'*w(:,1) + v'*Theta_mat_BS1_RIS_EU*w(:,1))^2 +...
                norm(ch_BS2_to_EU'*w(:,2)+ v'*Theta_mat_BS2_RIS_EU*w(:,2))^2 + noise_power) + (p(1)*(norm(ch_SU1_to_EU)^2 + ...
                norm(v'*theta_vec_SU1_RIS_EU)^2) + p(2)*(norm(ch_SU2_to_EU)^2 + norm(v'*theta_vec_SU2_RIS_EU)^2))/noise_power);
             
            seqobj(iIteration+1) = SU1_rate + SU2_rate;                
            seqobj(1:iIteration+1);

            if(abs(seqobj(iIteration+1)-seqobj(iIteration)) < tol)
% %                 keyboard;
                seqobj(iIteration+1:end)=[];
                v = v';
                status_cvx = 0;
                fprintf('PHASE-OPTI(CASE-3):%s\n', cvx_status);
                return;
            end            
        else
            fprintf('PHASE-OPTICASE-3):%s\n', cvx_status);
            status_cvx = 1;
            return;
        end
    end
end