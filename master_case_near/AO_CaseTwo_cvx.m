function [sumrate, SU1_rate, SU2_rate, EU_rate,pow_vec,v]=AO_CaseTwo_cvx(vec_theta, w, r_min_EU, PS1, PS2, PB1, PB2,...
            ch_BS1_to_EU, ch_BS2_to_EU, ch_RIS1_to_EU, ch_RIS2_to_EU, ch_BS1_to_RIS1, ch_BS2_to_RIS2, ch_SU1_to_EU, ...
            ch_SU2_to_EU, ch_BS1_to_SU1,  ch_BS2_to_SU2, ch_RIS1_to_SU1, ch_SU1_to_RIS1, ch_RIS2_to_SU2, ch_SU2_to_RIS2, noise_power, nIterations)
        
Theta_mat_BS1_RIS_SU1 = diag(ch_RIS1_to_SU1')*ch_BS1_to_RIS1;
Theta_mat_BS1_RIS_EU = diag(ch_RIS1_to_EU')*ch_BS1_to_RIS1;
theta_vec_SU1_RIS_EU = diag(ch_RIS1_to_EU')*ch_SU1_to_RIS1;

Theta_mat_BS2_RIS_SU2 = diag(ch_RIS2_to_SU2')*ch_BS2_to_RIS2;
Theta_mat_BS2_RIS_EU = diag(ch_RIS2_to_EU')*ch_BS2_to_RIS2;
theta_vec_SU2_RIS_EU = diag(ch_RIS2_to_EU')*ch_SU2_to_RIS2;

tol = 1e-3; 
 while(1)     
    seqobj = zeros(nIterations,1);
    for iIteration=1:nIterations % for AO algorithm 
        [w, pow_vec, cvx_status] = iterative_beamform_and_power_opti_case2_cvx(vec_theta, r_min_EU, PS1, PS2, PB1, PB2,...
                ch_BS1_to_EU, ch_BS2_to_EU, ch_RIS1_to_EU, ch_RIS2_to_EU, ch_BS1_to_RIS1, ch_BS2_to_RIS2, ch_SU1_to_EU, ...
                ch_SU2_to_EU, ch_BS1_to_SU1,  ch_BS2_to_SU2, ch_RIS1_to_SU1, ch_SU1_to_RIS1, ch_RIS2_to_SU2, ch_SU2_to_RIS2, noise_power, nIterations);                      
                   
        if (cvx_status ~= 0)
            sumrate = 0;
            SU1_rate = 0;
            SU2_rate = 0;
            EU_rate = 0;
            pow_vec = [0,0];
            v= zeros(200, 1);
            disp('Case 2: Failed');
            return;
        end
        [vec_theta, pow_vec, cvx_status] = iterative_phase_and_power_opti_case2_cvx(w, r_min_EU, PS1, PS2, PB1, PB2,...
                ch_BS1_to_EU, ch_BS2_to_EU, ch_RIS1_to_EU, ch_RIS2_to_EU, ch_BS1_to_RIS1, ch_BS2_to_RIS2, ch_SU1_to_EU, ...
                ch_SU2_to_EU, ch_BS1_to_SU1,  ch_BS2_to_SU2, ch_RIS1_to_SU1, ch_SU1_to_RIS1, ch_RIS2_to_SU2, ch_SU2_to_RIS2, noise_power, nIterations);         
        
        
                
        if (cvx_status ~= 0)
            sumrate = 0;
            SU1_rate = 0;
            SU2_rate = 0;
            EU_rate = 0;
            pow_vec = [0,0];
            v= zeros(200, 1);
            disp('Case 2: Failed');
            return;
        end
        norm([w(:,1); w(:,3)]);
        norm([w(:,2); w(:,4)]);
        % compute the Edge user's rate for each iteration
        v = vec_theta';
% compute the Edge user's rate for each iteration
        SU1_rate = log2(1 + norm(ch_BS1_to_SU1'*w(:,1) + v(:,1)'*Theta_mat_BS1_RIS_SU1*w(:,1))^2/noise_power);
        SU2_rate = log2(1 + norm(ch_BS2_to_SU2'*w(:,2) + v(:,2)'*Theta_mat_BS2_RIS_SU2*w(:,2))^2/noise_power);
        % compute the Edge user's rate for each iteration
        EU_rate_at_SU1 = log2(1 + norm(ch_BS1_to_SU1'*w(:,3) + v(:,1)'*Theta_mat_BS1_RIS_SU1*w(:,3))^2/(norm(ch_BS1_to_SU1'*w(:,1)+ v(:,1)'*Theta_mat_BS1_RIS_SU1*w(:,1))^2 + noise_power));
        EU_rate_at_SU2 = log2(1 + norm(ch_BS2_to_SU2'*w(:,4) + v(:,2)'*Theta_mat_BS2_RIS_SU2*w(:,4))^2/(norm(ch_BS2_to_SU2'*w(:,2)+ v(:,2)'*Theta_mat_BS2_RIS_SU2*w(:,2))^2 + noise_power));
        EU_rate = log2(1+ (norm(ch_BS1_to_EU'*w(:,3) + v(:,1)'*Theta_mat_BS1_RIS_EU*w(:,3))^2 + ...
            norm(ch_BS2_to_EU'*w(:,4)+ v(:,1)'*Theta_mat_BS2_RIS_EU*w(:,4))^2)/(norm(ch_BS1_to_EU'*w(:,1) + v(:,1)'*Theta_mat_BS1_RIS_EU*w(:,1))^2 +...
            norm(ch_BS2_to_EU'*w(:,2)+ v(:,2)'*Theta_mat_BS2_RIS_EU*w(:,2))^2 + noise_power) + (pow_vec(1)*(norm(ch_SU1_to_EU)^2 + ...
            norm(v(:,1)'*theta_vec_SU1_RIS_EU)^2) + pow_vec(2)*(norm(ch_SU2_to_EU)^2 + norm(v(:,2)'*theta_vec_SU2_RIS_EU)^2))/noise_power);
        
        seqobj(iIteration+1) = SU1_rate + SU2_rate;% + EU_rate;                
        seqobj(1:iIteration+1);

        if(abs(seqobj(iIteration+1)-seqobj(iIteration)) < tol)
%             keyboard;
            seqobj(iIteration+1:end)=[];
            beamformer_cvx = w;
            vec_theta_cvx = vec_theta;
            sumrate = seqobj(end);
            disp('Case 2: Successful');
            return;
        end 
    end
 end
        
