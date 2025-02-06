clear all

%% GENERAL
MU = 100;
LAM = 200;
A = 10;
N = 100;
TRIALS = 1e2;

%% CONFIG A1
% TAU = 1/sqrt(2*N);
% SEED = 1;
% FLAG_SIGMA_NORM = 0;
% SIGMA_LIST = linspace(0,1,21); SIGMA_LIST(1) = 1e-3; 

%% CONFIG A2
% TAU = 1/sqrt(2*N);
% SEED = 2;
% FLAG_SIGMA_NORM = 0;
% SIGMA_LIST = linspace(0,1,21); SIGMA_LIST(1) = 1e-3; 

%% CONFIG B1
% TAU = 1/sqrt(2*N);
% SEED = 1;
% FLAG_SIGMA_NORM = 1;
% SIGMA_LIST = linspace(0,50,21); SIGMA_LIST(1) = 1e-3; 

%% CONFIG B2
TAU = 1/sqrt(8*N);
SEED = 1;
FLAG_SIGMA_NORM = 1;
SIGMA_LIST = linspace(0,50,21); SIGMA_LIST(1) = 1e-3; 


%% Sphere
sph = Sphere;
[x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU,LAM,1,0), MU, LAM, N, 'full');
e_10 = e_vartheta_a_b(MU/LAM,1,0);
e_11 = e_vartheta_a_b(MU/LAM,1,1);

%% FIXED PARAMS (usually)
ALPHA = 2*pi;
N_G = 1;
SIGMA_STOP = nan;
F_STOP = nan;

%% FITNESS
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1); fit_name = 'Rastrigin';
ras_tools = RasTools;
Y_HAT = zeros(N, 1);

%% Y_0 INITIALIZATION
var = Variation;
[order_y_min, order_y_max] = ras_tools.get_all_extrema(ALPHA, A, N);

%% CREATE random Y with norm(Y)=R
POS_LIST = [sqrt(N)]; %sqrt(N)*[10,1,0.1,0.01]; %[100,10,1,0.1];
get_pos_var = @(y, n) var.get_Y_PPSN(POS_LIST, N, SEED);
Y_MAT = get_pos_var(POS_LIST, N);

%% CREATE Y=fixed
% SEED = 'nan';
% POS_LIST = [10];
% get_pos_var = @(y, n) var.get_Y_ones(POS_LIST, N);
% Y_MAT = get_pos_var(POS_LIST, N);

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
[psi_mean, psi_stderr] = deal(init_V1_V2, init_V1_V2);

%% RUN SIM
tic;
fprintf('START ...\n');

for p = 1:NUM_POS
    
    Y_0 = Y_MAT(:, p);
    fprintf('\n\t Y_0 = [%.2f, %.2f, ..., %.2f]', Y_0(1), Y_0(2), Y_0(end));
    
    for s = 1:NUM_SIG

        [phi_1_trial, phi_2_trial, phi_R_trial, phi_R2_trial, psi_trial] = ...
            deal(init_N_tr, init_N_tr, init_1_tr, init_1_tr, init_1_tr); 

        if FLAG_SIGMA_NORM==0
            SIGMA_0 = SIGMA_LIST(s);
            fprintf('\n\t\t Sigma = %.2f', SIGMA_0);
        else
            SIGMA_0 = SIGMA_LIST(s)*norm(Y_0)/N;
            fprintf('\n\t\t Sigma_norm=%.2f; (Sigma=%.2f)', SIGMA_LIST(s), SIGMA_0);
        end
        
        %ppm = ParforProgressbar(TRIALS);
        parfor t = 1:TRIALS
            
            F_0 = FIT(Y_0);
            R_0 = norm(Y_0 - Y_HAT);
            
            [y_opt, F_opt_vec, delta_r_vec, sigma_vec, gen] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, F_STOP, N_G, 0);
            
            phi_1_trial(:, t) = Y_0(:) - y_opt(:);
            phi_2_trial(:, t) = Y_0(:).^2 - y_opt(:).^2;
            phi_R_trial(1, t) = R_0 - norm(y_opt - Y_HAT);
            phi_R2_trial(1, t) = R_0^2 - norm(y_opt - Y_HAT)^2;
            psi_trial(1, t) = (sigma_vec(end) - SIGMA_0) / SIGMA_0;
            
            %ppm.increment();
        end % trials
        %delete(ppm);
        
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
        % Psi
        psi_mean(s, p) = mean(psi_trial);
        psi_stderr(s, p) = std(psi_trial, 0, 2)/sqrt(TRIALS);          
        
    end %sigma 

end %position
clear init_1_tr init_N_tr phi_1_R_trial psi_trial
clear phi_1_trial phi_2_trial phi_R_trial phi_R2_trial
fprintf('\nEND SIM\n'); 
toc;

%% Additional info
R = norm(Y_MAT(:,1));
fprintf('> Farthest minimum: %.2f \n', sqrt(N*order_y_min(end)^2));
title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES',', $\alpha$=$2\pi$', ', $A$=',num2str(A),', $N$=',num2str(N), ', $\tau$=', num2str(TAU),', $R$=',num2str(R), ', seed=',num2str(SEED)];
if FLAG_SIGMA_NORM==0; str_sig_label = '$\sigma$'; else, str_sig_label = '$\sigma^*$'; end

psi = Psi;
psi_res_y  = 0*SIGMA_LIST; dDQ=psi_res_y; DQ=psi_res_y; dEQ=psi_res_y;
psi_res_R  = 0*SIGMA_LIST;
psi_sph = 0*SIGMA_LIST; 
for i=1:length(SIGMA_LIST)
    if FLAG_SIGMA_NORM==0 
        sigma = SIGMA_LIST(i);
        signorm = SIGMA_LIST(i)*N/R;
    else
        sigma = SIGMA_LIST(i)*R/N;
        signorm = SIGMA_LIST(i);
    end
    psi_res_y(i) = psi_Y(A, ALPHA, TAU, sigma, Y_0, e_10, e_11);
    psi_res_R(i) = psi_R(A, ALPHA, TAU, sigma, R, N, e_10, e_11);
     %psi_sph(i) = psi.get_psi_sphere_signorm(MU, LAM, TAU, signorm);

end

figure; hold on; legend;
    title(title_str);
    xlabel(str_sig_label); ylabel('$\psi$');
    errorbar(SIGMA_LIST, psi_mean, psi_stderr, 'k.', 'DisplayName', 'Sim.');
    plot(SIGMA_LIST, psi_res_y, 'c-.', 'DisplayName', '$\psi(\mathbf{y})$');
    plot(SIGMA_LIST, psi_res_R, 'r--', 'DisplayName', '$\psi(R)$');
    %plot(SIGMA_LIST, psi_sph, 'k:', 'DisplayName', 'SPH');   
myfigstyle(gcf,6,5,8,8);
