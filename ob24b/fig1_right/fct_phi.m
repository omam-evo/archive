
function [phi_R_mean, phi_R_stderr] = fct_phi(MU, LAM, TAU, N, A, ALPHA, TRIALS, SIGMA_LIST)

%% PARAMS
% TRIALS = 1e4;
% A = 0;
% NUM_MU = length(MU_LIST);
% NUM_SIG = 11;

%% FIXED PARAMS (usually)
N_G = 1;
SIGMA_STOP = nan;
F_STOP = nan;

%% FITNESS
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1);

%% CREATE Y=const with R=const
rng(1);

%% Results
phi_R_mean = nan*SIGMA_LIST;
phi_R_stderr = nan*SIGMA_LIST;

%% RUN SIM
fprintf('START ...');
% figure; hold on;

    Y_0 = 1*ones(N,1);
    Y_HAT = zeros(N, 1);

    %
    NUM_POS = size(Y_0, 2);
    % Temp. variables for each variation
    init_1_tr = zeros(1, TRIALS);
    init_N_tr = zeros(N, TRIALS);

    fprintf('\n\t Y_0 = [%.2f, %.2f, ..., %.2f]', Y_0(1), Y_0(2), Y_0(end));

    FLAG_SIGMA_NORM = 1;

    for s = 1:length(SIGMA_LIST)

        [phi_1_trial, phi_2_trial, phi_R_trial, phi_R2_trial, q_gain_trial] = ...
            deal(init_N_tr, init_N_tr, init_1_tr, init_1_tr, init_1_tr); 

        if FLAG_SIGMA_NORM==0
            SIGMA_0 = SIGMA_LIST(s);
            fprintf('\n\t\t Sigma = %.2f', SIGMA_0);
        else
            SIGMA_0 = SIGMA_LIST(s)*norm(Y_0)/N;
            fprintf('\n\t\t Sigma_norm=%.2f; (Sigma=%.2f)', SIGMA_LIST(s), SIGMA_0);
        end
        
        for t = 1:TRIALS
            
            F_0 = FIT(Y_0);
            R_0 = norm(Y_0 - Y_HAT);
            
            [y_opt, F_opt_vec, ~, ~, gen] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, F_STOP, N_G, 0);

            % Phi_1, Phi_2, Phi_R, Q-Gain
            phi_1_trial(:, t) = Y_0(:) - y_opt(:);
            phi_2_trial(:, t) = Y_0(:).^2 - y_opt(:).^2;
            phi_R_trial(1, t) = R_0 - norm(y_opt);
            phi_R2_trial(1, t) = R_0^2 - norm(y_opt)^2;
            q_gain_trial(1, t)  = F_0 - F_opt_vec(gen+1);
            
        end % trials
        %delete(ppm);
        
        % Phi R
        phi_R_mean(s) = mean(phi_R_trial);
        phi_R_stderr(s) = std(phi_R_trial, 0, 2)/sqrt(TRIALS); 
        
    end %sigma

    phi_R_mean = phi_R_mean*N/R_0;
    phi_R_stderr = phi_R_stderr*N/R_0;

    fprintf('\nEND\n');
    
end


