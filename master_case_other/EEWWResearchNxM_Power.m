%------------------- (35)----(25-5-35)----------------
clc
clear all

tic;

% random seed
rng(5);

scenarios = 1;

M = 10; %number of antennas at BS
N = 200 ; %number of RIS elements
nIterations = 100; % No. of iteration needed for SCA algorithm

r_min_EU = 0.0001; 


%Base Station Power = Power at transmitter (Ptx)
pow_vec = [25]; %(2-14-22)

%Power variable for SUs
PS1_dBm = 25;%30 %35
PS1_dB = PS1_dBm-30;
PS1 = 10^(PS1_dB/10)  ;%PB1 to linear
PS2_dBm = 25;%30 %35
PS2_dB = PS2_dBm-30;
PS2 = 10^(PS2_dB/10)  ;%PB2 to linear

% noise power
noisePower_dB = -150; 
noise_power = 10^(noisePower_dB/10);

% Locations on X-Y coordinate (in meters)
BS1_loc = [-200, 0];
BS2_loc = [200, 0];

RIS1_loc = [-100,100];
RIS2_loc = [200,20];
RIS3_loc = [0,10];

edge_user_loc = [0,0];
cell_user_loc = [-120,20;100,30];

if 0
    figure
    plot(BS1_loc(:,1),BS1_loc(:,2),'r>');
    xlim([-210,210]);ylim([-110,110]);
    hold on;
    plot(BS2_loc(:,1),BS2_loc(:,2),'r>');
    plot(RIS1_loc(:,1),RIS1_loc(:,2),'bs');
    plot(RIS2_loc(:,1),RIS2_loc(:,2),'bs');
    plot(RIS3_loc(:,1),RIS3_loc(:,2),'bs');
    plot(cell_user_loc(:,1),cell_user_loc(:,2),'bo');
    plot(edge_user_loc(:,1),edge_user_loc(:,2),'ro');
end

% channel generation
eb=10;
eb2=1/(1+eb);
eb1=1-eb2;
eb1=sqrt(eb1);
eb2=sqrt(eb2);

BS_angle=rand(1,1);
RIS_angle=rand(1,1);

nRIS_rate_case1 = [];
nRIS_rate_case2 = [];
nRIS_rate_case3 = [];

disp('***Alternating Optimization Started***');
disp('--------------------------');
Sum_idx = 1:2 ;





for idx_npow = 1:length(pow_vec)
    fileID_case1_vs_pow = fopen("Sumrate_case1_vs_pow_" + pow_vec(idx_npow) + "dBm.txt",'w');

    fileID_case2_vs_pow = fopen("Sumrate_case2_vs_pow_" + pow_vec(idx_npow) + "dBm.txt",'w');
    
    fileID_case3_vs_pow = fopen("Sumrate_case3_vs_pow_" + pow_vec(idx_npow) + "dBm.txt",'w');
    
    fileID_p_case1         =fopen("P_case1_"+ pow_vec(idx_npow) + "dBm.txt",'w');
    fileID_v_case1         =fopen("V_case1_"+ pow_vec(idx_npow) + "dBm.txt",'w');

    fileID_p_case2         =fopen("P_case2_"+ pow_vec(idx_npow) + "dBm.txt",'w');
    fileID_v_case2         =fopen("V_case2_"+ pow_vec(idx_npow) + "dBm.txt",'w');

    fileID_p_case3       =fopen("P_case3_"+ pow_vec(idx_npow) + "dBm.txt",'w');
    fileID_v_case3       =fopen("V_case3_"+ pow_vec(idx_npow) + "dBm.txt",'w');

    PB1_dBm = pow_vec(idx_npow);%power of BS1   ---------------------------------------------------------------------
    PB1_dB = PB1_dBm-30;
    PB1 = 10^(PB1_dB/10)  ;%PB1 to linear
    
    PB2_dBm = pow_vec(idx_npow);%power of BS2   ----------------------------------------------------------------------
    PB2_dB = PB2_dBm-30;
    PB2 = 10^(PB2_dB/10)  ;%PB2 to linear

    
    total_rate_case1 = [];
    total_rate_case2 = [];
    total_rate_case3 = [];
    disp(['Number of RIS elements are:' num2str(N)]); 
    
    
    for idx=1:length(Sum_idx)
    %channel generation
    % ch_BS1_to_EU1, ch_BS_to_EU2, ch_RIS_to_EU1, ch_RIS_to_EU2, ch_BS1_to_RIS, ch_SU1_to_EU1, ch_SU2_to_EU2, ch_BS1_to_SU1,  ch_BS2_to_SU2
    
    % NLoS link: BS1 to EU
    dist_BS1_to_EU = sqrt((BS1_loc(1) - edge_user_loc(1))^2 + (BS1_loc(2) - edge_user_loc(2))^2);
    path_BS1_to_EU = 10^(-path_NLOS(dist_BS1_to_EU)/10);
    ch_BS1_to_EU1 = sqrt(path_BS1_to_EU).*(1/sqrt(2))*(rand(M,1)+rand(M,1)*1i);%channel from BS1 to EU
    % NLoS link: BS2 to EU
    dist_BS2_to_EU = sqrt((BS2_loc(1) - edge_user_loc(1))^2 + (BS2_loc(2) - edge_user_loc(2))^2);
    path_BS2_to_EU = 10^(-path_NLOS(dist_BS2_to_EU)/10);
    ch_BS2_to_EU1 = sqrt(path_BS2_to_EU).*(1/sqrt(2))*(rand(M,1)+rand(M,1)*1i);%channel from BS2 to EU
    % NLoS link: BS1 to SU1
    dist_BS1_to_SU1 = sqrt((BS1_loc(1) - cell_user_loc(1,1))^2 + (BS1_loc(2) - cell_user_loc(1,2))^2);
    path_BS1_to_SU1 = 10^(-path_NLOS(dist_BS1_to_SU1)/10);
    ch_BS1_to_SU1 = sqrt(path_BS1_to_SU1).*(1/sqrt(2))*(rand(M,1)+rand(M,1)*1i);%channel from BS1 TO SU1
    % NLoS link: BS2 to SU2
    dist_BS2_to_SU2 = sqrt((BS2_loc(1) - cell_user_loc(2,1))^2 + (BS2_loc(2) - cell_user_loc(2,2))^2);
    path_BS2_to_SU2 = 10^(-path_NLOS(dist_BS2_to_SU2)/10);
    ch_BS2_to_SU2 = sqrt(path_BS2_to_SU2).*(1/sqrt(2))*(rand(M,1)+rand(M,1)*1i);%channel from BS2 TO SU2
    % NLoS link: SU1 to EU
    dist_SU1_to_EU = sqrt((cell_user_loc(1,1) - edge_user_loc(1))^2 + (cell_user_loc(1,2) - edge_user_loc(2))^2);
    path_SU1_to_EU = 10^(-path_NLOS(dist_SU1_to_EU)/10);
    ch_SU1_to_EU1 = sqrt(path_SU1_to_EU).*(1/sqrt(2))*(rand(1,1)+rand(1,1)*1i);% channel from SU1 to EU
    % NLoS link: SU2 to EU
    dist_SU2_to_EU = sqrt((cell_user_loc(2,1) - edge_user_loc(1))^2 + (cell_user_loc(2,2) - edge_user_loc(2))^2);
    path_SU2_to_EU = 10^(-path_NLOS(dist_SU2_to_EU)/10);
    ch_SU2_to_EU1 = sqrt(path_SU2_to_EU).*(1/sqrt(2))*(rand(1,1)+rand(1,1)*1i);%channel from SU2 to EU
    % LOS link: RIS1 to EU
    dist_RIS1_to_EU = sqrt((RIS1_loc(1) - edge_user_loc(1))^2 + (RIS1_loc(2) - edge_user_loc(2))^2);
    path_RIS1_to_EU = 10^(-path_LOS(dist_RIS1_to_EU)/10);
    ch_RIS1_to_EU_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_RIS1_to_EU1 = sqrt(path_RIS1_to_EU)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_RIS1_to_EU_NLoS);%channel from RIS to EU
    % LOS link: RIS2 to EU
    dist_RIS2_to_EU = sqrt((RIS2_loc(1) - edge_user_loc(1))^2 + (RIS2_loc(2) - edge_user_loc(2))^2);
    path_RIS2_to_EU = 10^(-path_LOS(dist_RIS2_to_EU)/10);
    ch_RIS2_to_EU_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_RIS2_to_EU2 = sqrt(path_RIS2_to_EU)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_RIS2_to_EU_NLoS);%channel from BS2 RIS to EU 
    % LOS link: RIS3 to EU
    dist_RIS3_to_EU = sqrt((RIS3_loc(1) - edge_user_loc(1))^2 + (RIS3_loc(2) - edge_user_loc(2))^2);
    path_RIS3_to_EU = 10^(-path_LOS(dist_RIS3_to_EU)/10);
    ch_RIS3_to_EU_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_RIS3_to_EU = sqrt(path_RIS3_to_EU)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_RIS3_to_EU_NLoS);
    % LOS link: RIS1 to SU1
    dist_RIS1_to_SU1 = sqrt((RIS1_loc(1) - cell_user_loc(1,1))^2 + (RIS1_loc(2) - cell_user_loc(1,2))^2);
    path_RIS1_to_SU1 = 10^(-path_LOS(dist_RIS1_to_SU1)/10);
    ch_RIS1_to_SU1_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_RIS1_to_SU1 = sqrt(path_RIS1_to_SU1)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_RIS1_to_SU1_NLoS);%channel from BS1 RIS to EU 
    % LOS link: RIS2 to SU2
    dist_RIS2_to_SU2 = sqrt((RIS2_loc(1) - cell_user_loc(2,1))^2 + (RIS2_loc(2) - cell_user_loc(2,2))^2);
    path_RIS2_to_SU2 = 10^(-path_LOS(dist_RIS2_to_SU2)/10);
    ch_RIS2_to_SU2_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_RIS2_to_SU2 = sqrt(path_RIS2_to_SU2)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_RIS2_to_SU2_NLoS);
    % LOS link: RIS3 to SU1
    dist_RIS3_to_SU1 = sqrt((RIS3_loc(1) - cell_user_loc(1,1))^2 + (RIS3_loc(2) - cell_user_loc(1,2))^2);
    path_RIS3_to_SU1 = 10^(-path_LOS(dist_RIS3_to_SU1)/10);
    ch_RIS3_to_SU1_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_RIS3_to_SU1 = sqrt(path_RIS3_to_SU1)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_RIS3_to_SU1_NLoS);
    % LOS link: RIS3 to SU2
    dist_RIS3_to_SU2 = sqrt((RIS3_loc(1) - cell_user_loc(2,1))^2 + (RIS3_loc(2) - cell_user_loc(2,2))^2);
    path_RIS3_to_SU2 = 10^(-path_LOS(dist_RIS3_to_SU2)/10);
    ch_RIS3_to_SU2_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_RIS3_to_SU2 = sqrt(path_RIS3_to_SU2)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_RIS3_to_SU2_NLoS);
    % LOS link: SU1 to RIS1
    ch_SU1_to_RIS1_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_SU1_to_RIS1 = sqrt(path_RIS1_to_SU1)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_SU1_to_RIS1_NLoS);
    % LOS link: SU2 to RIS2
    ch_SU2_to_RIS2_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_SU2_to_RIS2 = sqrt(path_RIS2_to_SU2)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_SU2_to_RIS2_NLoS);
    % LOS link: SU1 to RIS3
    ch_SU1_to_RIS3_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_SU1_to_RIS3 = sqrt(path_RIS3_to_SU1)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_SU1_to_RIS3_NLoS);
    % LOS link: SU2 to RIS3
    ch_SU2_to_RIS3_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_SU2_to_RIS3 = sqrt(path_RIS3_to_SU2)*(eb1.*ULA_fun(RIS_angle ,N)+eb2.*ch_SU2_to_RIS3_NLoS);
    % LOS link: BS1 to RIS1
    dist_BS1_to_RIS1 = sqrt((BS1_loc(1) - RIS1_loc(1))^2 + (BS1_loc(2) - RIS1_loc(2))^2);
    path_BS1_to_RIS1 = 10^(-path_LOS(dist_BS1_to_RIS1)/10);
    ch_BS1_to_RIS1_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_BS1_to_RIS1 = sqrt(path_BS1_to_RIS1)*(eb1.*ULA_fun(RIS_angle ,N)*ULA_fun(BS_angle ,M)'+eb2.*ch_BS1_to_RIS1_NLoS);%channel from BS1 to RIS1 
    % LOS link: BS2 to RIS2
    dist_BS2_to_RIS2 = sqrt((BS2_loc(1) - RIS2_loc(1))^2 + (BS2_loc(2) - RIS2_loc(2))^2);
    path_BS2_to_RIS2 = 10^(-path_LOS(dist_BS2_to_RIS2)/10);
    ch_BS2_to_RIS2_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_BS2_to_RIS2 = sqrt(path_BS2_to_RIS2)*(eb1.*ULA_fun(RIS_angle ,N)*ULA_fun(BS_angle ,M)'+eb2.*ch_BS2_to_RIS2_NLoS);%channel from BS2 to RIS2 ****************
    % LOS link: BS1 to RIS3
    dist_BS1_to_RIS3 = sqrt((BS1_loc(1) - RIS3_loc(1))^2 + (BS1_loc(2) - RIS3_loc(2))^2);
    path_BS1_to_RIS3 = 10^(-path_LOS(dist_BS1_to_RIS3)/10);
    ch_BS1_to_RIS3_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_BS1_to_RIS3 = sqrt(path_BS1_to_RIS3)*(eb1.*ULA_fun(RIS_angle ,N)*ULA_fun(BS_angle ,M)'+eb2.*ch_BS1_to_RIS3_NLoS);
    % LOS link: BS2 to RIS3
    dist_BS2_to_RIS3 = sqrt((BS2_loc(1) - RIS3_loc(1))^2 + (BS2_loc(2) - RIS3_loc(2))^2);
    path_BS2_to_RIS3 = 10^(-path_LOS(dist_BS2_to_RIS3)/10);
    ch_BS2_to_RIS3_NLoS = (1/sqrt(2))*(randn(N, 1)+j*randn(N, 1));
    ch_BS2_to_RIS3 = sqrt(path_BS2_to_RIS3)*(eb1.*ULA_fun(RIS_angle ,N)*ULA_fun(BS_angle ,M)'+eb2.*ch_BS2_to_RIS3_NLoS);
    
    w_int = random_unit_vector(M,4);

%         case 1
    vec_theta_case1 = rand(1,N)*pi/2;
    phase_vec_case1 = exp(1j*vec_theta_case1.*2.*pi);
%     [sumrate_case1, SU1_rate_case1, SU2_rate_case1, EU_rate_case1,p,v]=AO_CaseOne_cvx(phase_vec_case1, w_int, r_min_EU, PS1, PS2, PB1, PB2,...
%     ch_BS1_to_EU1, ch_BS2_to_EU1, ch_RIS1_to_EU1, ch_RIS2_to_EU2, ch_BS1_to_RIS1, ch_SU1_to_EU1, ...
%     ch_SU2_to_EU1, ch_BS1_to_SU1,  ch_BS2_to_SU2, ch_RIS1_to_SU1, ch_SU1_to_RIS1, noise_power, nIterations);
% sumSEsystem(idx_npow,idx,:) = [pow_vec(idx_npow), idx, SU1_rate_case1, SU2_rate_case1, EU_rate_case1];
%             fprintf(fileID_case1_vs_pow,'%7d %10d %10.2f %10.4f %10.2f\n',sumSEsystem(idx_npow,idx,:));
%             
%        % case 2
%     vec_theta_case2 = rand(2,N)*pi/2;
%     phase_vec_case2 = exp(1j*vec_theta_case2.*2.*pi);
%     [sumrate_case2, SU1_rate_case2, SU2_rate_case2, EU_rate_case2,p2,v2]=AO_CaseTwo_cvx(phase_vec_case2, w_int, r_min_EU, PS1, PS2, PB1, PB2,...
%     ch_BS1_to_EU1, ch_BS2_to_EU1, ch_RIS1_to_EU1, ch_RIS2_to_EU2, ch_BS1_to_RIS1, ch_BS2_to_RIS2, ch_SU1_to_EU1, ...
%     ch_SU2_to_EU1, ch_BS1_to_SU1,  ch_BS2_to_SU2, ch_RIS1_to_SU1, ch_SU1_to_RIS1, ch_RIS2_to_SU2, ch_SU2_to_RIS2, noise_power, nIterations);
% sumSEsystem(idx_npow,idx,:) = [pow_vec(idx_npow), idx, SU1_rate_case2, SU2_rate_case2, EU_rate_case2];
%             fprintf(fileID_case2_vs_pow,'%7d %10d %10.2f %10.4f %10.2f\n',sumSEsystem(idx_npow,idx,:));

%         case 3
    [sumrate_case3, SU1_rate_case3, SU2_rate_case3, EU_rate_case3,p3,v3]=AO_CaseThree_cvx(phase_vec_case1, w_int, r_min_EU, PS1, PS2, PB1, PB2,...
    ch_BS1_to_EU1, ch_BS2_to_EU1, ch_RIS3_to_EU, ch_BS1_to_RIS3, ch_BS2_to_RIS3, ch_SU1_to_EU1, ...
    ch_SU2_to_EU1, ch_BS1_to_SU1,  ch_BS2_to_SU2, ch_RIS3_to_SU1, ch_SU1_to_RIS3, ch_RIS3_to_SU2, ch_SU2_to_RIS3, noise_power, nIterations);
sumSEsystem(idx_npow,idx,:) = [pow_vec(idx_npow), idx, SU1_rate_case3, SU2_rate_case3, EU_rate_case3];
            fprintf(fileID_case3_vs_pow,'%7d %10d %10.2f %10.4f %10.2f\n',sumSEsystem(idx_npow,idx,:));

% total_rate_case1(idx) = SU1_rate_case1 + SU2_rate_case1 + EU_rate_case1;
% total_rate_case2(idx) = SU1_rate_case2 + SU2_rate_case2 + EU_rate_case2;
total_rate_case3(idx) = SU1_rate_case3 + SU2_rate_case3 + EU_rate_case3;
% y = [idx, total_rate_case1(idx), SU1_rate_case1, SU2_rate_case1, EU_rate_case1];
% y2 = [idx, total_rate_case2(idx), SU1_rate_case2, SU2_rate_case2, EU_rate_case2];
y3 = [idx, total_rate_case3(idx), SU1_rate_case3, SU2_rate_case3, EU_rate_case3];
disp('**************************');
disp('    Ch. #  | Sumrate | SU1 rate | SU2 rate | EU rate')
% disp(y)
% disp(y2)
disp(y3)
disp('**************************');

%     total_rate;

    end
%     nRIS_rate_case1(idx_npow) = mean(total_rate_case1)
%     SU_case1(idx_npow)        = mean(SU1_rate_case1 + SU2_rate_case1)
%     EU_case1(idx_npow)        = mean (EU_rate_case1)
%     nRIS_rate_case2(idx_npow) = mean(total_rate_case2)
%     SU_case2(idx_npow)        = mean(SU1_rate_case2 + SU2_rate_case2)
%     EU_case2(idx_npow)        = mean (EU_rate_case2)
    nRIS_rate_case3(idx_npow) = mean(total_rate_case3)
    SU_case3(idx_npow)        = mean(SU1_rate_case3 + SU2_rate_case3)
    EU_case3(idx_npow)        = mean (EU_rate_case3)
% fprintf(fileID_p_case1,'%7.3d  \n',p);
% fprintf(fileID_v_case1,'%7.3d  \n',v);
% 
% fprintf(fileID_p_case2,'%7.3d  \n',p2);
% fprintf(fileID_v_case2,'%7.3d  \n',v2);

fprintf(fileID_p_case3,'%7.3d  \n',p3);
fprintf(fileID_v_case3,'%7.3d  \n',v3);
end
% sumrate_case1 = mean(nRIS_rate_case1);
% sumrate_case2 = mean(nRIS_rate_case2);
sumrate_case3 = mean(nRIS_rate_case3);
toc
