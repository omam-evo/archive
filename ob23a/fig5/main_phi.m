clear

%% CONFIG 1: const Y
% MODE = 'constY';

%% CONFIG 2: const R
MODE = 'constR';

%% Params
MU = 100;
LAM = 200;
A = 1;
N = 100;
ALPHA = 2*pi;
N_G = 1;
TRIALS = 1e4;

%% Position
SEED = 1;
rng(SEED);
R = 7;
Y_MAT = R/sqrt(N)*ones(N,1);
POS_LIST = [R]; 

%% SIGMA Regular
% FLAG_SIGMA_NORM = 0;
% SIGMA_LIST = [0:0.1:1]'; SIGMA_LIST(1) = 1e-5; 
%% SIGMA Norm
FLAG_SIGMA_NORM = 1;
SIGMA_LIST = linspace(0,50,26)'; SIGMA_LIST(1) = 1e-5; 

%% 1 Generation experiments
TAU = 0;
SIGMA_STOP = nan;
F_STOP = nan;

%% FITNESS
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1);
Y_HAT = zeros(N, 1);

%% DATA INITIALIZATION
NUM_SIG = length(SIGMA_LIST);
NUM_POS = size(Y_MAT, 2);
% Temp. variables for each variation
init_1_tr = zeros(1, TRIALS);
init_N_tr = zeros(N, TRIALS);
% Persistent variable: avg. over trials
init_V1_V2 = zeros(NUM_SIG, NUM_POS);
init_N_V1_V2 = zeros(N, NUM_SIG, NUM_POS);
% Init
[phi_1_mean, phi_1_stderr] = deal(init_N_V1_V2, init_N_V1_V2);
[phi_2_mean, phi_2_stderr] = deal(init_N_V1_V2, init_N_V1_V2);
[phi_R_mean, phi_R_stderr] = deal(init_V1_V2, init_V1_V2);
[phi_R2_mean, phi_R2_stderr] = deal(init_V1_V2, init_V1_V2);
[q_gain_mean]  = deal(init_V1_V2);

%% RUN SIM
tic;
fprintf('START ...\n');

for p = 1:NUM_POS
    
    % Y_0 is overwritten below for mode constR
    Y_0 = Y_MAT(:, p);
    if strcmpi(MODE, 'constY')
        fprintf('\n\t Y_0 = [%.2f, %.2f, ..., %.2f]', Y_0(1), Y_0(2), Y_0(end));
    elseif strcmpi(MODE, 'constR')
        fprintf('\n\t R_0 = %.2f', POS_LIST(p));
    end
    
    for s = 1:NUM_SIG

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
            
            % change Y_0 value with constR
            if strcmpi(MODE, 'constR')
                v = randn(N,1);
                Y_0 = v/norm(v)*POS_LIST(p); %override const value Y_0
            end
            
            F_0 = FIT(Y_0);
            R_0 = norm(Y_0 - Y_HAT);
            
            [y_opt, F_opt_vec, delta_r_vec, sigma_vec, gen] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, F_STOP, N_G, 0);
            
            % Phi_1, Phi_2, Phi_R, Q-Gain
            phi_1_trial(:, t) = Y_0(:) - y_opt(:);
            phi_2_trial(:, t) = Y_0(:).^2 - y_opt(:).^2;
            phi_R_trial(1, t) = R_0 - norm(y_opt - Y_HAT);
            phi_R2_trial(1, t) = R_0^2 - norm(y_opt - Y_HAT)^2;
            q_gain_trial(1, t)  = F_0 - F_opt_vec(gen+1);
            
        end % trials
        
        % Phi 1
        phi_1_mean(:, s, p) = mean(phi_1_trial, 2);
        phi_1_stderr(:, s, p) = std(phi_1_trial, 0, 2)/sqrt(TRIALS); 
        % Phi 2
        phi_2_mean(:, s, p) = mean(phi_2_trial, 2);
        phi_2_stderr(:, s, p) = std(phi_2_trial, 0, 2)/sqrt(TRIALS);  
        % Phi R
        phi_R_mean(s, p) = mean(phi_R_trial);
        phi_R_stderr(s, p) = std(phi_R2_trial, 0, 2)/sqrt(TRIALS); 
        phi_R2_mean(s, p) = mean(phi_R2_trial);
        phi_R2_stderr(s, p) = std(phi_R2_trial, 0, 2)/sqrt(TRIALS);  
        % Quality Gain
        q_gain_mean(s, p) = mean(q_gain_trial, 2);
        
    end %sigma 

end %position
clear init_1_tr init_N_tr phi_1_R_trial q_gain_trial
clear phi_1_trial phi_2_trial phi_R_trial phi_R2_trial
fprintf('\nEND SIM\n'); 
toc;

save('wsp.mat');

           
%% RUN APPROXIMATION
FLAG_NORM_PHI = 1;
SIGMA_LIST_highRes = linspace(SIGMA_LIST(1),SIGMA_LIST(end),101);
NUM_SIG = length(SIGMA_LIST_highRes);
phi_R2_B2 = zeros(NUM_SIG, NUM_POS);
phi_R2 = zeros(NUM_SIG, NUM_POS);

for p=1:NUM_POS
    Y_0 = Y_MAT(:, p);
    %Y_i = Y_0(COMP_ID); Y_jni = Y_0(COMP_ID ~= 1:N);
    R = norm(Y_0);
    N = length(Y_0);
    if FLAG_NORM_PHI == 1
        phi_R2_reg2norm = N/(2*R^2);
    else
        phi_R2_reg2norm = 1;
    end
    phi_R2_mean(:,p) = phi_R2_mean(:,p)*phi_R2_reg2norm;
    phi_R2_stderr(:,p) = phi_R2_stderr(:,p)*phi_R2_reg2norm;

    for s=1:NUM_SIG

        if FLAG_SIGMA_NORM==0, sigma = SIGMA_LIST_highRes(s); else, sigma = SIGMA_LIST_highRes(s)*norm(Y_0)/N; end

        %% Phi 2
        phi_R2_B2(s, p) = phi_R2_reg2norm * sum(phi_2(MU, LAM, A, ALPHA, sigma, Y_0, 0, 0),1);
        phi_R2(s, p) = phi_R2_reg2norm * phi_R(MU, LAM, N, A, ALPHA, sigma, R);

    end
    
    %% Plot
    figure; hold on;
        if FLAG_SIGMA_NORM==0; str_sig_label = '$\sigma$'; else, str_sig_label = '$\sigma^*$'; end
        xlabel(str_sig_label); 
        ylabel('$\varphi_R^{\mathrm{II},*}$');

        errorbar(SIGMA_LIST, phi_R2_mean(:, p), phi_R2_stderr(:, p), 'k.', 'DisplayName', 'SIM');
        if strcmpi(MODE, 'constY')
            plot(SIGMA_LIST_highRes, phi_R2_B2(:, p), 'b-.', 'DisplayName', 'B2');
        elseif strcmpi(MODE, 'constR')
            plot(SIGMA_LIST_highRes, phi_R2(:, p), 'r--', 'DisplayName', 'R');
        end
    xlim([0,SIGMA_LIST(end)]);
    ylim([-2,8]);
    myfigstyle(gcf,7.5,4,9,9);

end
