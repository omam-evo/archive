

%% CONFIG 1; iteration with small pop. and N
% MU = 10;
% LAM = 40;
% A = 1;
% N = 20;
% TRIALS = 1e6;
% POS_LIST = [sqrt(N)*10,sqrt(N)*0.1]; %[sqrt(N)*10,sqrt(N)*1,sqrt(N)*0.1,sqrt(N)*0.01];
% FLAG_SIGMA_NORM = 1;
% SIGMA_LIST = linspace(0,20,21)'; SIGMA_LIST(1) = 1e-5; 

%% CONFIG 2; iteration with larger pop. and N
MU = 100;
LAM = 200;
A = 1;
N = 100;
TRIALS = 1e3;
POS_LIST = [sqrt(N)*10,sqrt(N)*0.1];
FLAG_SIGMA_NORM = 1;
SIGMA_LIST = linspace(0,60,21)'; SIGMA_LIST(1) = 1e-5; 

file = 'main_phi';
SAVE = 0;
sph = Sphere;
[signorm_opt, signorm_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});

%% SIGMA Regular
% FLAG_SIGMA_NORM = 0;
% SIGMA_LIST = [0:0.05:1]'; SIGMA_LIST(1) = 1e-5; 
%% SIGMA Norm
% sigmax = 10+ceil(signorm_2ndZero/10)*10;
% FLAG_SIGMA_NORM = 1;
% SIGMA_LIST = linspace(0,sigmax,21)'; SIGMA_LIST(1) = 1e-5; 

%% FIXED PARAMS (usually)
ALPHA = 2*pi;
N_G = 1;
TAU = 0;
SIGMA_STOP = nan;
F_STOP = nan;

%% FITNESS
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1); fit_name = 'Rastrigin';
ras_tools = RasTools;
Y_HAT = zeros(N, 1);

%% Y_0 INITIALIZATION
var = Variation;
[order_y_min, order_y_max] = ras_tools.get_all_extrema(ALPHA, A, N);

%% CREATE Y BY R
POS_LABEL = '$R$';
SEED = 1;
% POS_LIST = [sqrt(N)*1]; %[sqrt(N)*10,sqrt(N)*1,sqrt(N)*0.1,sqrt(N)*0.01];
get_pos_var = @(y, n) var.get_Y_PPSN(POS_LIST, N, SEED);
Y_MAT = get_pos_var(POS_LIST, N);

%% CREATE Y
% POS_LABEL = '$Y_0$';
% SEED = 1;
% POS_LIST = [1.25]; %[100,10,1,0.1];
% get_pos_var = @(y, n) var.get_Y_ones(POS_LIST, N, SEED);
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
[q_gain_mean]  = deal(init_V1_V2);

%% RUN SIM
tic;
fprintf('START ...\n');

for p = 1:NUM_POS
    
    Y_0 = Y_MAT(:, p);
    fprintf('\n\t Y_0 = [%.2f, %.2f, ..., %.2f]', Y_0(1), Y_0(2), Y_0(end));
    
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
        
        %ppm = ParforProgressbar(TRIALS);
        parfor t = 1:TRIALS

            F_0 = FIT(Y_0);
            R_0 = norm(Y_0 - Y_HAT);
            
            %% Default RNG
            [y_opt, F_opt_vec, delta_r_vec, sigma_vec, gen] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, F_STOP, N_G, 0);
            
            % Phi_1, Phi_2, Phi_R, Q-Gain
            phi_1_trial(:, t) = Y_0(:) - y_opt(:);
            phi_2_trial(:, t) = Y_0(:).^2 - y_opt(:).^2;
            phi_R_trial(1, t) = R_0 - norm(y_opt - Y_HAT);
            phi_R2_trial(1, t) = R_0^2 - norm(y_opt - Y_HAT)^2;
            q_gain_trial(1, t)  = F_0 - F_opt_vec(gen+1);
            
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
        % Quality Gain
        q_gain_mean(s, p) = mean(q_gain_trial, 2);
        
    end %sigma 

end %position
clear init_1_tr init_N_tr phi_1_R_trial q_gain_trial
clear phi_1_trial phi_2_trial phi_R_trial phi_R2_trial
fprintf('\nEND SIM\n'); 
toc;

%% Additional info
fprintf('> Farthest minimum: %.2f \n', sqrt(N*order_y_min(end)^2));
fprintf('> Sphere signorm 2nd zero: %.2f \n', signorm_2ndZero);

if SAVE == 1
    savestr = ['M',num2str(MU),'_L',num2str(LAM),'_V',num2str(MU/LAM),'_N',num2str(N),'_A',num2str(A)];
    save(fullfile(path_folder,['wsp_',savestr,'_',datestr(now,'yyyy-mm-dd_HH-MM-SS'),'.mat'])); 
end 
% end %main ===== DOE =====

%% LOAD DATA and execute the evaluation part below

%% INIT
sph = Sphere;
ras_tools = RasTools;
ep = EvalPhi;
FLAG_NORM_PHI = 1;

%% Phi 1
COMP_ID = 2; getpath = '';
[~, phi_1_A3] = ep.get_phi1_i_plot(ras_tools, MU, LAM, A, ALPHA, SIGMA_LIST, Y_MAT, ...
    phi_1_mean, phi_1_stderr, FLAG_SIGMA_NORM, COMP_ID, getpath, {'SIM','A3'}); %'ki', 'di'

myfigstyle(gcf, 18, 4, 9, 9);

%% Phi 2
sph = Sphere;
ras_tools = RasTools;
ep = EvalPhi;
FLAG_NORM_PHI = 1;
COMP_ID = 2; getpath = '';
[phi_2_num, phi_2_B1, phi_2_B2, phi_2_L1, phi_2_L2] = ep.get_phi2_i(ras_tools, MU, LAM, A, ALPHA, SIGMA_LIST, Y_MAT, phi_2_mean, phi_2_stderr, FLAG_SIGMA_NORM, COMP_ID, ...
    getpath, {'SIM','B2','L1','PLOT'});
myfigstyle(gcf, 18, 5, 9, 9);

%% Error measure
% er = Error;
% sigma_get = 10;
% id = find(SIGMA_LIST==sigma_get);
% if ~isempty(id)
%     [error_A3] = er.difference_at_location(phi_1_mean, phi_1_A3);
%     figure; hold on; xlabel('Component $i$'); ylabel('$\varphi_{i, approx}-\varphi_{i, meas}$');
%         title(['$\sigma=',num2str(SIGMA_LIST(id)), '$']); legend;
%         plot(1:N, error_A3(:,id), 'b.', 'DisplayName', 'A3');
%         plot([1,N], [mean(error_A3(:,id)),mean(error_A3(:,id))], 'b--', 'DisplayName', 'mean A3')
%         myfigstyle(gcf, 16, 10,10,10);
% else
%     warning('NOT FOUND: SIGMA_LIST==sigma_get')
% end
% myfigstyle(gcf, 9, 5,9,9);
