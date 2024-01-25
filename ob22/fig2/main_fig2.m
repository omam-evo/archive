clear
% addpath('../parforProgress')

%% VARIATION
R = 100; %values R=100,10,1,0.1

%% SETTINGS
HOLD_ON = 0; % legend('AutoUpdate','off');
SUBPLOT = 0;
SIZE_X = 15;
SIZE_Y = 10;

%% PARAMS
MU = 150;
LAM = 300;
% A = 10; %see sigma
N = 100;
TRIALS = 1e5;

%% FIXED PARAMS (usually)
ALPHA = 2*pi;
N_G = 1;
TAU = 0;
SIGMA_STOP = nan;
F_STOP = nan;

%% SIGMA CONF_1 ==> ylim([0,0.03]); two error config = {0.1,1}
% PLOT = 1;
% A = 10;
% R = 10;
% FLAG_SIGMA_NORM = 0;
% SIGMA_LIST = linspace(0, 1, 21)';
% SIGMA_LIST(1) = 1e-5; 
% COMP_ID = 2; % COMP=2 => y_i=1.19; COMP=12 => y_i=0.8

%% SIGMA_NORM CONF_2 ==> ylim([0,])
PLOT = 2;
A = 1;
FLAG_SIGMA_NORM = 1;
SIGMA_LIST = [0:2.5:40]';
SIGMA_LIST(1) = 1e-5; 
COMP_ID = 2;

%% FITNESS
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1); fit_name = 'Rastrigin';
ras_tools = RasTools;
Y_HAT = zeros(N, 1);

%% Y_0 INITIALIZATION
var = Variation;
[order_y_min, order_y_max] = ras_tools.get_all_extrema(ALPHA, A, N);
%% ONES
% Y_0 = 1.25*ones(N,1);
% Y_0(2:end) = 0; %Y_0(floor(N/2):end);
%% Random
% rand('seed', 10); Y_0 = 3*rand(N,1); Y_0(1) = 1.71;
% rand('seed', 3); Y_0 = 3*rand(N,1); Y_0(1) = 1.25;
%% MIN/MAX
% Y_0 = order_y_min(2)*ones(N,1);
% Y_0 = order_y_max(2)*ones(N,1);

%% CREATE Y BY ZEROS
% POS_LABEL = '$Y_0(1)$, $Y_0$(2:end)=0';
% POS_START = 10; POS_END = 10; NUM_POS = 1;
% POS_LIST = var.get_position_list(POS_START, POS_END, NUM_POS);
% get_pos_var = @(a, b) var.get_Y_zerosExcept1(a, b, nan);
% Y_MAT = get_pos_var(POS_LIST, N);

%% CREATE Y BY ONES
% POS_LABEL = '$Y_0$*ones$(N,1)$';
% POS_START = 1.75; POS_END = 1.75; NUM_POS = 1;
% POS_LIST = var.get_position_list(POS_START, POS_END, NUM_POS);
% get_pos_var = @(y, n) var.get_Y_ones(y, n, nan);
% Y_MAT = get_pos_var(POS_LIST, N);

%% CREATE Y BY R
POS_LABEL = '$R$';
POS_START = R; POS_END = R; NUM_POS = 1; SEED = 1; %sqrt(N*10^2)
POS_LIST = var.get_position_list(POS_START, POS_END, NUM_POS);
get_pos_var = @(y, n) var.get_Y_from_R(y, n, SEED);
Y_MAT = get_pos_var(POS_LIST, N);

%% DATA INITIALIZATION
NUM_SIG = length(SIGMA_LIST);
NUM_POS = size(Y_MAT, 2);
% Temp. variables for each variation
init_1_tr = zeros(1, TRIALS);
init_N_tr = zeros(N, TRIALS);
% Persistent variable: avg. over trials
init_1_V1_V2 = zeros(1, NUM_SIG, NUM_POS);
init_N_V1_V2 = zeros(N, NUM_SIG, NUM_POS);
% Init
[phi_1_mean, phi_1_stderr] = deal(init_N_V1_V2, init_N_V1_V2);
[phi_2_mean, phi_2_stderr] = deal(init_N_V1_V2, init_N_V1_V2);
[phi_1_R_mean] = deal(init_1_V1_V2);
[q_gain_mean]  = deal(init_1_V1_V2);

%% RUN SIM
tic;
fprintf('START ...\n');

for p = 1:NUM_POS
    
    Y_0 = Y_MAT(:, p);
    fprintf('\n\t Y_0 = [%.2f, %.2f, ..., %.2f]', Y_0(1), Y_0(2), Y_0(end));
    
    for s = 1:NUM_SIG

        [phi_1_trial, phi_2_trial, phi_1_R_trial, q_gain_trial] = ...
            deal(init_N_tr, init_N_tr, init_1_tr, init_1_tr); 

        if FLAG_SIGMA_NORM==0
            SIGMA_0 = SIGMA_LIST(s);
            fprintf('\n\t\t Sigma = %.2f', SIGMA_0);
        else
            SIGMA_0 = SIGMA_LIST(s)*norm(Y_0)/N;
            fprintf('\n\t\t Sigma_norm=%.2f; (Sigma=%.2f)', SIGMA_LIST(s), SIGMA_0);
        end
        
%         ppm = ParforProgressbar(TRIALS);
        parfor t = 1:TRIALS

            rstream = RandStream('threefry4x64_20');
            rstream.Substream = t;

            F_0 = FIT(Y_0);
            R_0 = norm(Y_0 - Y_HAT);

            [y_opt, F_opt_vec, delta_r_vec, sigma_vec, gen] = muComLam_sigSA_randStream(rstream, FIT, N, MU, LAM,Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, F_STOP, N_G, 0);
            %[y_opt, F_opt_vec, delta_r_vec, sigma_vec, gen] = muComLam_sSA(FIT, N, MU, LAM,Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, F_STOP, N_G, 0);

            % Phi_1, Phi_2, Phi_R, Q-Gain
            phi_1_trial(:, t) = Y_0(:) - y_opt(:);
            phi_2_trial(:, t) = Y_0(:).^2 - y_opt(:).^2;
            phi_1_R_trial(1, t) = R_0 - norm(y_opt - Y_HAT);
            q_gain_trial(1, t)  = F_0 - F_opt_vec(gen+1);
            
%             ppm.increment();
        end % trials
%         delete(ppm);
        
        % Phi 1
        phi_1_mean(:, s, p) = mean(phi_1_trial, 2);
        phi_1_stderr(:, s, p) = std(phi_1_trial, 0, 2)/sqrt(TRIALS); 
        % Phi 2
        phi_2_mean(:, s, p) = mean(phi_2_trial, 2);
        phi_2_stderr(:, s, p) = std(phi_2_trial, 0, 2)/sqrt(TRIALS);  
        % Phi R
        phi_1_R_mean(1, s, p) = mean(phi_1_R_trial, 2);
        % Quality Gain
        q_gain_mean(1, s, p) = mean(q_gain_trial, 2);
        
    end %sigma 

end %position
clear init_1_tr init_N_tr phi_1_R_trial q_gain_trial
clear phi_1_trial phi_2_trial
fprintf('\nEND SIM\n'); 
toc;

%% INIT SIG V2
FACTOR_APPROX = 1;
NUM_SIG_V2 = NUM_SIG*FACTOR_APPROX;
SIGMA_LIST_V2 = linspace(SIGMA_LIST(1), SIGMA_LIST(end), NUM_SIG_V2);

%% INIT POS V2
NUM_POS_V2 = NUM_POS*FACTOR_APPROX;
POS_LIST_V2 = var.get_position_list(POS_START, POS_END, NUM_POS_V2);
Y_MAT_V2 = get_pos_var(POS_LIST_V2, N);

%% INIT PHI
Phi = Phi_v2;
C_MU_LAM = e_mu_lam_a_b(MU, LAM, 1, 0); % only needed for phi_1 via c_mu_lam
phi_1_A3 = zeros(N, NUM_SIG_V2, NUM_POS_V2);
phi_1_C_k = zeros(N, NUM_SIG_V2, NUM_POS_V2);
phi_1_C_f = zeros(N, NUM_SIG_V2, NUM_POS_V2);

%% INIT PHI_2
% e_11 = e_mu_lam_a_b(MU, LAM, 1, 1);       % only needed for phi_2
% e_20 = e_mu_lam_a_b(MU, LAM, 2, 0);       % only needed for phi_2
% phi_2_B1 = zeros(N, NUM_SIG_V2, NUM_POS_V2);
% phi_2 = zeros(N, NUM_SIG_V2, NUM_POS_V2);
        
%% RUN APPROXIMATION
for p=1:NUM_POS_V2
    Y_0 = Y_MAT_V2(:, p);

    for s=1:NUM_SIG_V2

        if FLAG_SIGMA_NORM==0, sig = SIGMA_LIST_V2(s); else, sig = SIGMA_LIST_V2(s)*norm(Y_0)/N; end
        
        [phi, ~, ~] = Phi.get_phi_1_viaC_Di(ras_tools, C_MU_LAM, A, ALPHA, sig, Y_0, 0);
        phi_1_C_k(:, s, p) = phi;
        [phi, ~, ~] = Phi.get_phi_1_viaC_Di(ras_tools, C_MU_LAM, A, ALPHA, sig, Y_0, 1);
        phi_1_C_f(:, s, p) = phi;
        [phi, ~, ~] = Phi.get_phi_1(ras_tools, MU, LAM, A, ALPHA, sig, Y_0, 1);
        phi_1_A3(:, s, p) = phi;

    end
end

%% Error measure
if PLOT == 1
    er = Error;
    % [error, error_id, phi_error] = er.mean_squared(phi_2_mean, phi_2);
    sigma_get = 1;
    id = find(SIGMA_LIST==sigma_get);
    [error_C_k] = er.difference_at_location(phi_1_mean, phi_1_C_k);
    [error_C_f] = er.difference_at_location(phi_1_mean, phi_1_C_f);
    figure; hold on; xlabel('Component $i$'); ylabel('$\varphi_{i, approx}-\varphi_{i, meas}$');
        title(['$\sigma=',num2str(SIGMA_LIST(id)), '$']); legend;
        % plot(1:N, error_phi_2, 'go')
        % yline(mean(error_phi_2), 'g--')
        plot(1:N, error_C_k(:,id), 'b.', 'DisplayName', '$k_i$')
        %yline(mean(error_C_k), 'r--', 'DisplayName', 'mean')
        plot(1:N, error_C_f(:,id), 'r.', 'DisplayName', '$f_i''$')
        %yline(mean(error_C_f), 'b--', 'DisplayName', 'mean')
        myfigstyle(gcf, 16, 10,10,10);
end

%% Figure
POS_ID = 1;
Y_0 = Y_MAT(:, POS_ID);

if FLAG_SIGMA_NORM==0; str_sig_label = '$\sigma$'; else, str_sig_label = '$\sigma^*$'; end

if HOLD_ON == 0
    Y_0 = Y_MAT_V2(:, POS_ID);
    Y_i = Y_0(COMP_ID); 
    Y_inj = Y_0(COMP_ID ~= 1:N);
    Y_str = []; %,', Y_j=[', num2str(Y_inj(1),3),', ',num2str(Y_inj(2),3),', ...,',num2str(Y_inj(end),3), ']']; 
    title_str = ['$A=',num2str(A),',~R=',num2str(R),',~y_i=', num2str(Y_i,3), '$'];
    if SUBPLOT == 1
        figure('units', 'centimeters', 'position', [0, 0, 32, 16]);
        sgtitle(title_str, 'interpreter', 'latex'); figuresize(32, 16, 'cm');
    else
        figure(); title(title_str, 'interpreter', 'latex')
        figuresize(SIZE_X, SIZE_Y, 'cm');
    end
else
    legend("AutoUpdate", 'off');
end

%% Q-Gain and Phi_R
if SUBPLOT == 1
subplot(2, 3, 1); hold on;%subplot(2,2,1); 

    %% Q-Gain
    xlabel(str_sig_label,'fontsize',14,'interpreter','latex');

    %% Phi_R
    yyaxis left;
        ax = gca; ax.YColor = 'k';
        ylabel('$\varphi_R = R^{(g)}-R^{(g+1)}$','fontsize',14,'interpreter','latex');
        plot(SIGMA_LIST, phi_1_R_mean(1, :, POS_ID), 'ko:', 'DisplayName', 'Exp.');
        yline(0, 'k--');
        % plot(sigma_list, mean(phi_2_all), 'ro:', 'DisplayName', 'Exp.');

    yyaxis right; 
        ylabel('$Q = f^{(g)}-f^{(g+1)}$','fontsize',14,'interpreter','latex'); % [pos.=imprv.]
        ax = gca; ax.YColor = 'b';
        plot(SIGMA_LIST, q_gain_mean(1, :, POS_ID), 'bo:', 'DisplayName', 'Exp.');
        yline(0, 'b--');  

end %folding

%% Components
if SUBPLOT == 1
    subplot(2, 3, [2,3,5,6]);      
end
    hold on;
    xlabel(str_sig_label,'fontsize',18,'interpreter','latex'); ylabel('$\varphi_i$','fontsize',18,'interpreter','latex');

        %% Phi experimental
        errorbar(SIGMA_LIST, phi_1_mean(COMP_ID, :, POS_ID), phi_1_stderr(COMP_ID, :, POS_ID), 'k.', 'DisplayName', 'Sim.', 'LineWidth', 1, 'MarkerSize', 15);
        % errorbar(SIGMA_LIST, phi_2_mean(COMP_ID, :, POS_ID), phi_2_stderr(COMP_ID, :, POS_ID), 'k.', 'DisplayName', '$\varphi_i^{\mathrm{II}}$', 'LineWidth', 1, 'MarkerSize', 15);

        %% Numerical exact
        % plot(SIGMA_LIST, phi_i_NUM, 'k:', 'DisplayName', '$\varphi_i$: NUM', 'LineWidth', 2);

        %% Analytical approx.
        %plot(SIGMA_LIST_V2, phi_1_A3(COMP_ID, :, POS_ID), 'b--', 'DisplayName', '$\varphi_i(k_i + \exp[..]d_i)$', 'LineWidth', 2);
        plot(SIGMA_LIST_V2, phi_1_C_k(COMP_ID, :, POS_ID), 'b-.', 'DisplayName', '$\varphi_i(k_i)$', 'LineWidth', 2);
        plot(SIGMA_LIST_V2, phi_1_C_f(COMP_ID, :, POS_ID), 'r--', 'DisplayName', '$\varphi_i(f_i'')$', 'LineWidth', 2);
        
        %plot(SIGMA_LIST_V2, phi_2(COMP_ID, :, POS_ID), 'r-.', 'DisplayName', '$\varphi_i^{\mathrm{II}} = 2y_i\varphi_i^{(\mathrm{A3})} - \sigma^2/\mu$', 'LineWidth', 2);
        %plot(SIGMA_LIST_V2, phi_2_B1(COMP_ID, :, POS_ID), 'g-.', 'DisplayName', '$\varphi_i^{\mathrm{II}} = 2y_i\varphi_i^{(\mathrm{A3})} - \sigma^2/\mu - [...]$', 'LineWidth', 2);

        %% Update
        legend('location', 'best', 'interpreter', 'latex');
        
    %% Rastrigin
    if SUBPLOT==1
        subplot(2, 3, 4);hold on; %subplot(2,2,3); 
        Y_i = Y_0(COMP_ID);
        xlabel('$y_i$','fontsize',12,'interpreter','latex'); ylabel('$f(y_i)$','fontsize',12,'interpreter','latex');
        xrange = linspace(Y_i-1.5, Y_i+1.5, 100);
        plot(xrange, FIT(xrange), '-', 'Color', [0.3,0.3,0.3]);
        plot(Y_i, FIT(Y_i), 'ko', 'MarkerSize', 10);
        xline(Y_i, '-', 'Color', [0.5,0.5,0.5]);
        syms x; T = ras_tools.get_taylor(ALPHA, A, 1, Y_i);
        fplot(T(x),[Y_i-0.2, Y_i+0.2], 'r--', 'DisplayName', 'Taylor: linear');
    end % folding

myfigstyle(gcf, 16, 10, 10, 10);