function [v, p, status_cvx] = iterative_phase_and_power_opti_cvx(w, r_min, PS1, PS2, PB1, PB2,...
            ch_BS1_to_EU, ch_BS2_to_EU, ch_RIS_to_EU, ch_RIS_to_EU2, ch_BS1_to_RIS, ch_SU1_to_EU, ...
            ch_SU2_to_EU, ch_BS1_to_SU1,  ch_BS2_to_SU2, ch_RIS_to_SU1, ch_SU1_to_RIS, noise_power, nIterations)
%         noise_power = 1;

[N,M] = size(ch_BS1_to_RIS);
Theta_mat_BS1_RIS_SU1 = diag(ch_RIS_to_SU1')*ch_BS1_to_RIS;
Theta_mat_BS1_RIS_EU = diag(ch_RIS_to_EU')*ch_BS1_to_RIS;
theta_vec_SU1_RIS_EU = diag(ch_RIS_to_EU')*ch_SU1_to_RIS;

r_min_EU = 2^r_min -1;
% t_next = 0.005*ones(1,1);
p_next = 0.00005*ones(1,1);
eta_next = 0.00005*ones(3,1);
pt_next = 0.000005*ones(4,1);
qt_next = 0.000005*ones(4,1);
x_next = 0.00005*ones(N,1);
y_next = 0.00005*ones(N,1);

scale_1 = 1e-2/noise_power;
tol = 1e-2;     
% cvx_solver Mosek
 cvx_solver sedumi
%cvx_solver 
cvx_precision high
% cvx_solver_settings('MSK_IPAR_INFEAS_REPORT_AUTO','MSK_ON','MSK_IPAR_INFEAS_REPORT_LEVEL',1)
while(1)
    sumrate = zeros(nIterations,1);
    seqobj = zeros(nIterations,1);
    for iIteration=1:nIterations % for SCA algorithm
        cvx_begin quiet  
%         variables   t(2) z eta(3) p(2) u b(2) pt(3) qt(3);
        variables   x(N,1) y(N,1) z(N,1) t b eta(2) pt(4) qt(4) p(2) u theta(N,1);
        variable    v(N,1) complex;

        maximize( t + sum( x_next.^2 + y_next.^2 + 2.*x_next.*(x-x_next) + 2*y_next.*(y-y_next)));
         
        subject to
        for ii=1:N
            x(ii) == real(v(ii));
            y(ii) == imag(v(ii));
        end
        % UE1 Eq. 19d       
        norm([sqrt(noise_power); ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1)]) <= sqrt(eta(1)); % Eq. 30d        
        
        % Edge UE Eq. 19g
        norm([sqrt(noise_power); ch_BS2_to_EU'*w(:,2); ch_BS1_to_EU'*w(:,1) + v'*Theta_mat_BS1_RIS_EU*w(:,1)]) <= sqrt(eta(2)); %Eq. 30g
        
        % Eq. 19e
        pt(4) == real(v'*theta_vec_SU1_RIS_EU);
        qt(4) == imag(v'*theta_vec_SU1_RIS_EU);
        scale_1*(p(1)*norm(ch_SU1_to_EU)^2 + p(1)*(pt_next(4)^2 + qt_next(4)^2) + 2*p_next(1)*pt_next(4)*(pt(4) - pt_next(4))...
           +2*p_next(1)*qt_next(4)*(qt(4) - qt_next(4)) + p(2)*norm(ch_SU2_to_EU)^2 + noise_power*u) >= scale_1*noise_power*r_min_EU;
       
       % Eq. 19h
        for i=1:N
                norm([x(i); y(i)]) <= 1;
        end
        % Eq. 19i & 19j 
        p(1) <= PS1; % Eq. 30i
        p(2) <= PS2; % Eq. 30j
        
        
        % Eq. 19b
        pt(1) == real(ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1));
        qt(1) == imag(ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1));
        scale_1*(pt_next(1)*(2*pt(1)-pt_next(1)) + qt_next(1)*(2*qt(1) - qt_next(1))) >= scale_1*(t-1)*noise_power;
%         2*(pt_next(1)/t_next)*(pt(1)-pt_next(1)) + 2*(qt_next(1)/t_next)*(qt(1) - qt_next(1))  + ((pt_next(1)^2 + qt_next(1)^2)/(t_next^2))*(1-(t-t_next)/t_next)>= noise_power;
        
       
%         % Eq. 19c
        pt(2) == real(ch_BS1_to_SU1'*w(:,3) + v'*Theta_mat_BS1_RIS_SU1*w(:,3));
        qt(2) == imag(ch_BS1_to_SU1'*w(:,3) + v'*Theta_mat_BS1_RIS_SU1*w(:,3));
        pt_next(2)*(2*pt(2)-pt_next(2)) + qt_next(2)*(2*qt(2) - qt_next(2)) >= eta(1)*r_min_EU ;         

     
        % Eq. 19f
        pt(3) == real(ch_BS1_to_EU'*w(:,3) + v'*Theta_mat_BS1_RIS_EU*w(:,3));
        qt(3) == imag(ch_BS1_to_EU'*w(:,3) + v'*Theta_mat_BS1_RIS_EU*w(:,3));
        2*pt_next(3)*(pt(3)-pt_next(3))/eta_next(2) + 2*qt_next(3)*(qt(3) - qt_next(3))/eta_next(2) + ((pt_next(3)^2 + qt_next(3)^2)/eta_next(2))*(1-(eta(2) - eta_next(2))/eta_next(2)) + (norm(ch_BS2_to_EU'*w(:,4))^2/eta_next(2))*(2-eta(2)/eta_next(2)) >= u  ;
        

        p(1) >= 0;
        p(2) >= 0;
        t >=  1 ;
        u >= 0;
        eta(1) >= 0;
        eta(2) >= 0;
  
        cvx_end;
        if(strfind(cvx_status,'Solved'))
%             t_next = double(t);
            p_next = p;
            eta_next = double(eta);
            pt_next = double(pt);
            qt_next = double(qt);
            x_next = x;
            y_next = y;

            
            norm(w(:,1))^2 + norm(w(:,3))^2;
            norm(w(:,2))^2 + norm(w(:,4))^2;
            % compute the Edge user's rate for each iteration
            SU1_rate = log2(1 + norm(ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1))^2/noise_power);
            SU2_rate = log2(1+norm(ch_BS2_to_SU2'*w(:,2))^2/noise_power);
%             % compute the Edge user's rate for each iteration
            EU_rate_at_SU1 = log2(1 + norm(ch_BS1_to_SU1'*w(:,3)...
                + v'*Theta_mat_BS1_RIS_SU1*w(:,3))^2/(norm(ch_BS1_to_SU1'*w(:,1)...
                + v'*Theta_mat_BS1_RIS_SU1*w(:,1))^2 + noise_power));
            EU_rate_at_SU2 = log2(1 + norm(ch_BS2_to_SU2'*w(:,4))^2/(norm(ch_BS2_to_SU2'*w(:,2))^2 + noise_power));
            EU_rate = log2(1+ (norm(ch_BS1_to_EU'*w(:,3) + v'*Theta_mat_BS1_RIS_EU*w(:,3))^2 + ...
                norm(ch_BS2_to_EU'*w(:,4))^2)/(norm(ch_BS1_to_EU'*w(:,1) + v'*Theta_mat_BS1_RIS_EU*w(:,1))^2 +...
                norm(ch_BS2_to_EU'*w(:,2))^2 + noise_power) + (p(1)*(norm(ch_SU1_to_EU)^2 + ...
                norm(v'*theta_vec_SU1_RIS_EU)^2) + p(2)*norm(ch_SU2_to_EU)^2)/noise_power);
             
            seqobj(iIteration+1) = SU1_rate + SU2_rate;                
            seqobj(1:iIteration+1);

            if(abs(seqobj(iIteration+1)-seqobj(iIteration)) < tol)
% %                 keyboard;
                seqobj(iIteration+1:end)=[];
                v = v';
                status_cvx = 0;
                fprintf('PHASE-OPTI(CASE-1):%s\n', cvx_status');
                return;
            end            
        else
            fprintf('PHASE-OPTI(CASE-1):%s\n', cvx_status');
            status_cvx = 1;
            return;
        end
    end
end