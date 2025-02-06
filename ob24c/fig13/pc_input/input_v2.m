function [pc_input] = input_v2(pc_input, N) %input_slow_rescale

    % Timescale
    TS = max(0, N^(1/2));
    cumulation = 1/TS;
    pc_input.wait = 0;
    
    % Sigma
    pc_input.sigma_rescale = 0;
    
    % alpha_mu
    pc_input.fac_inc = 1.05;
    pc_input.fac_dec = 1/pc_input.fac_inc;

    % additional settings
    pc_input.pid_method = 'pid_v1'; % Pop inc/dec: 'pid_v1', 'pid_v2', 'pid_psa'
    pc_input.damp = 0; % used by 'pid_v2'
    pc_input.annotate = 0;
    
    %% Methods
    % 'PCCMSA', 'APOP'
    pc_input.L = ceil(TS); 
    % 'APOP'
    pc_input.apop_thresh = 0.2;
    % 'PCCMSA'
    pc_input.alpha_sign = 0.05; 
    % 'PSA'
    pc_input.alpha_psa = 1.4; %1.4;
    pc_input.beta_psa = cumulation; %0.4;

    %% Test
    % % 'ROSAES'
    % pc_input.c_f = cumulation; %'ROSAES'
    % pc_input.robust_thresh = 0; %'ROSAES'
    % % 'RANK'
    % pc_input.alpha_rank = 0.05;
    % % 'LENGTHZ'
    % pc_input.alpha_lengthz = 1-0.9545; % 2sigma
    % pc_input.c_lengthz = cumulation;
    
    %% Work in progress
    % 'INFNOISE'
    % pc_input.c_infnoise = cumulation; 
    %p c_input.stddev_infnoise = @(MU,N) sqrt(2/(MU*N)*(1-1/MU));
    % pc_input.ev_infnoise = @(MU) (1-1/MU);
    % pc_input.alpha_infnoise = 0.05;
    % % 'QGAIN'
    % pc_input.qgain_thresh = 0; 
    % pc_input.c_qgain = 1;

end