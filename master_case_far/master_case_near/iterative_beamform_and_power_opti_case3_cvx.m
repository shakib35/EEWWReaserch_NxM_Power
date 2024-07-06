function [w, p, status_cvx] = iterative_beamform_and_power_opti_case3_cvx(vec_theta, r_min, PS1, PS2, PB1, PB2,...
            ch_BS1_to_EU, ch_BS2_to_EU, ch_RIS3_to_EU, ch_BS1_to_RIS3, ch_BS2_to_RIS3, ch_SU1_to_EU, ...
                ch_SU2_to_EU, ch_BS1_to_SU1,  ch_BS2_to_SU2, ch_RIS3_to_SU1, ch_SU1_to_RIS3, ch_RIS3_to_SU2, ch_SU2_to_RIS3, noise_power, nIterations)
%         noise_power = 1;
v = vec_theta';
[N,M] = size(ch_BS1_to_RIS3);
Theta_mat_BS1_RIS_SU1 = diag(ch_RIS3_to_SU1')*ch_BS1_to_RIS3;
Theta_mat_BS1_RIS_EU = diag(ch_RIS3_to_EU')*ch_BS1_to_RIS3;
theta_vec_SU1_RIS_EU = diag(ch_RIS3_to_EU')*ch_SU1_to_RIS3;

Theta_mat_BS2_RIS_SU2 = diag(ch_RIS3_to_SU2')*ch_BS2_to_RIS3;
Theta_mat_BS2_RIS_EU = diag(ch_RIS3_to_EU')*ch_BS2_to_RIS3;
theta_vec_SU2_RIS_EU = diag(ch_RIS3_to_EU')*ch_SU2_to_RIS3;

EU_rate = 0;

scale1= 1000;
r_min_EU = 2^r_min -1;
t_next = 0.005*ones(2,1);
p_next = 0.0005*ones(2,1);
eta_next = 0.0005*ones(3,1);
pt_next = 0.005*ones(8,1);
qt_next = 0.005*ones(8,1);
u_next = 0.005*ones(1,1);
tol = 1e-2;     
cvx_solver Mosek
cvx_precision high

while(1)
    sumrate = zeros(nIterations,1);
    seqobj = zeros(nIterations,1);
    for iIteration=1:nIterations
        cvx_begin quiet  
        variables   t(2) z eta(3) p(2) u b(2) pt(8) qt(8);
        variable    w(M,4) complex; 
        
        maximize(z);
         
        subject to
        % Objective function
        norm([b(1)-b(2); 2*z]) <= b(1)+b(2);       
        
        % UE1  Eq. 43f      
        norm([sqrt(noise_power); ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1)]) <= sqrt(eta(1)); 
         
        % UE2 Eq. 43g
        norm([sqrt(noise_power); ch_BS2_to_SU2'*w(:,2) + v'*Theta_mat_BS2_RIS_SU2*w(:,2)]) <= sqrt(eta(2)); 
        
        % Edge UE Eq. 43j
        norm([sqrt(noise_power); ch_BS2_to_EU'*w(:,2)+ v'*Theta_mat_BS2_RIS_EU*w(:,2); ch_BS1_to_EU'*w(:,1) + v'*Theta_mat_BS1_RIS_EU*w(:,1)]) <= sqrt(eta(3));
        

        % Eq. 43h
        pt(6) >= real(v'*theta_vec_SU1_RIS_EU);
        qt(6) >= imag(v'*theta_vec_SU1_RIS_EU);
        pt(8) >= real(v'*theta_vec_SU2_RIS_EU);
        qt(8) >= imag(v'*theta_vec_SU2_RIS_EU);
        p(1)*(pt_next(6)^2 + qt_next(6)^2) + 2*p_next(1)*pt_next(6)*(pt(6) - pt_next(6)) +2*p_next(1)*qt_next(6)*(qt(6) - qt_next(6))...
       +p(2)*(pt_next(8)^2 + qt_next(8)^2) + 2*p_next(2)*pt_next(8)*(pt(8) - pt_next(8)) +2*p_next(2)*qt_next(8)*(qt(8) - qt_next(8))...
       >= noise_power*r_min_EU - noise_power*u - p(1)*norm(ch_SU1_to_EU)^2  - p(2)*norm(ch_SU2_to_EU)^2 ;
         
        % Eq. 43k & 43l
        norm([w(:,1); w(:,3)]) <= sqrt(PB1);
        norm([w(:,2); w(:,4)]) <= sqrt(PB2);

        % Eq. 43m & 43n
        p(1) <= PS1;
        p(2) <= PS2;
               
        t(1) >= b(1)-1;
        t(2) >= b(2)-1;
        % Eq. 43b
        pt(1) == (real(ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1)));
        qt(1) == (imag(ch_BS1_to_SU1'*w(:,1) + v'*Theta_mat_BS1_RIS_SU1*w(:,1)));
        2*pt_next(1)*(pt(1)-pt_next(1))/t_next(1) + 2*qt_next(1)*(qt(1) - qt_next(1))/t_next(1) + ((pt_next(1)^2 + qt_next(1)^2)/t_next(1))*(1-(t(1) - t_next(1))/t_next(1))>= noise_power;
%         pt_next(1)*(2*pt(1)-pt_next(1)) + qt_next(1)*(2*qt(1) - qt_next(1)) >= t(1)*noise_power;
        
        % Eq. 43c 
        pt(2) == (real(ch_BS2_to_SU2'*w(:,2) + v'*Theta_mat_BS2_RIS_SU2*w(:,2)));
        qt(2) == (imag(ch_BS2_to_SU2'*w(:,2) + v'*Theta_mat_BS2_RIS_SU2*w(:,2)));
        2*pt_next(2)*(pt(2)-pt_next(2))/t_next(2) + 2*qt_next(2)*(qt(2) - qt_next(2))/t_next(2) + ((pt_next(2)^2 + qt_next(2)^2)/t_next(2))*(1-(t(2) - t_next(2))/t_next(2))>= noise_power;
%         pt_next(2)*(2*pt(2)-pt_next(2)) + qt_next(2)*(2*qt(2) - qt_next(2)) >= t(2)*noise_power;
        
        % Eq. 43d
        pt(3) == real(ch_BS1_to_SU1'*w(:,3) + v'*Theta_mat_BS1_RIS_SU1*w(:,3));
        qt(3) == imag(ch_BS1_to_SU1'*w(:,3) + v'*Theta_mat_BS1_RIS_SU1*w(:,3));
        2*pt_next(3)*(pt(3)-pt_next(3))/eta_next(1) + 2*qt_next(3)*(qt(3) - qt_next(3))/eta_next(1) + ((pt_next(3)^2 + qt_next(3)^2)/eta_next(1))*(1-(eta(1) - eta_next(1))/eta_next(1))>= r_min_EU;
        
        % Eq. 43e
        pt(4) == real(ch_BS2_to_SU2'*w(:,4)+ v'*Theta_mat_BS2_RIS_SU2*w(:,4));
        qt(4) == imag(ch_BS2_to_SU2'*w(:,4)+ v'*Theta_mat_BS2_RIS_SU2*w(:,4));
        2*pt_next(4)*(pt(4)-pt_next(4))/eta_next(2) + 2*qt_next(4)*(qt(4) - qt_next(4))/eta_next(2) + ((pt_next(4)^2 + qt_next(4)^2)/eta_next(2))*(1-(eta(2) - eta_next(2))/eta_next(2))>= r_min_EU;
         
        % Eq. 43i
        pt(5) == real(ch_BS1_to_EU'*w(:,3) + v'*Theta_mat_BS1_RIS_EU*w(:,3));
        qt(5) == imag(ch_BS1_to_EU'*w(:,3) + v'*Theta_mat_BS1_RIS_EU*w(:,3));
        pt(7) == real(ch_BS2_to_EU'*w(:,4) + v'*Theta_mat_BS2_RIS_EU*w(:,4));
        qt(7) == imag(ch_BS2_to_EU'*w(:,4) + v'*Theta_mat_BS2_RIS_EU*w(:,4));
        2*pt_next(5)*(pt(5)-pt_next(5))/eta_next(3) + 2*qt_next(5)*(qt(5) - qt_next(5))/eta_next(3) + ((pt_next(5)^2 + qt_next(5)^2)/eta_next(3))*(1-(eta(3) - eta_next(3))/eta_next(3)) +...
        2*pt_next(7)*(pt(7)-pt_next(7))/eta_next(3) + 2*qt_next(7)*(qt(7) - qt_next(7))/eta_next(3) + ((pt_next(7)^2 + qt_next(7)^2)/eta_next(3))*(1-(eta(3) - eta_next(3))/eta_next(3)) >= u;
        
        p(1) >= 0;
        p(2) >= 0;
        b(1) >= 1;
        b(2) >= 1;
        z >= 0;
        t(1) >= 0;
        t(2) >= 0;
        u >= 0;
        eta >= 0;
        cvx_end;
        w;
        if(strfind(cvx_status,'Solved'))
            
            t_next = double(t);
            p_next = p;
            eta_next = double(eta);
            pt_next = double(pt);
            qt_next = double(qt);
            u_next = double(u);
            
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
%                 keyboard;
                seqobj(iIteration+1:end)=[];
                status_cvx = 0;
                fprintf('BEAMFORM-OPTI(CASE-3):%s\n', cvx_status');
                return;
            end
        else
            fprintf('BEAMFORM-OPTICASE-3):%s\n', cvx_status');
            status_cvx = 1;
            return;
        end
    end
end  