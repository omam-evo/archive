clear
file = 'main_dyn';

%% SETTINGS
SAVE = 0;

%% CONFIG 1
MU = 10;
LAM = 40;
A = 1;
N = 20;
ALPHA = 2*pi;
TRIALS = 100;
N_G = 1000;
SIGMA_NORM = 7;

%% CONFIG 2
% MU = 100;
% LAM = 200;
% A = 1;
% N = 100;
% ALPHA = 2*pi;
% TRIALS = 100;
% N_G = 1000;
% SIGMA_NORM = 30;

C_MU_LAM = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
sph = Sphere; 
% SIGMA_NORM = 1.0*sph.signorm_opt_sph(C_MU_LAM, MU, LAM, N, {'medium'});
SIGMA_STOP = 1e-10; %1e-10;
R_STOP = 1e-8;
F_STOP = nan;
tau_str = '1/sqrt(2*N)';

%% Plot options
PLOT_AVG = 1;
PLOT_ALL_R = 1;

%% Position fix
% Y_0 = 10*ones(N, 1); R_0 = norm(Y_0);

%% Position random
a = 20;   % R > 7*sqrt(N) for A=10 necessary
R_0 = a*sqrt(N);
SEED = 1; rng(SEED); v = randn(N,1); 
Y_0 = v/norm(v)*R_0;

%% Set SIGMA & TAU
% if MU == 1, TAU = 1/sqrt(N_DIM); else, TAU = 1/sqrt(2*N_DIM); end
SIGMA_0 = MU*e_vartheta_a_b(MU/LAM, 1, 0)*R_0/N;
if isnan(SIGMA_NORM)
    TAU = eval(tau_str);  %smaller tau for better convergence
    sig_ctrl_str = ['$\tau$=', tau_str, ''];
else
    sig_ctrl_str = ['$\sigma^*$=', num2str(SIGMA_NORM)];
end

%% Fitness
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1); fit_name = 'Rastrigin';
Y_HAT = zeros(N, 1);

%% Get all extrema
ras_tools = RasTools;
[order_y_min, order_y_max] = ras_tools.get_all_extrema(ALPHA, A, N);

%% Initialize
init_vec = zeros(N_G+1, 1); % first value contains 0-th generation
init_trial = zeros(TRIALS, 1);
[sigma_avg, sigma_norm_avg] = deal(init_vec, init_vec);
[r_avg, f_avg] = deal(init_vec, init_vec);
[flag_success, num_gen, num_feval] = deal(init_trial, init_trial, init_trial);
r_g_t = zeros(N_G+1, TRIALS); 

%% Main Loop 
tic
for i = 1:TRIALS

    fprintf('  t = %i \n', i)

    if isnan(SIGMA_NORM)
        [y_opt, F_opt_g, r_g, sigma_g, gen] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, N_G, 0);
    else
        [y_opt, F_opt_g, r_g, sigma_g, gen] = muComLam_sigNorm(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_NORM, SIGMA_STOP, R_STOP, N_G, 0);
    end
    
    %% Assertion failes if N_G is reached before F_STOP or SIGMA_STOP => increase max. value of N_G
    %assert(sigma_g(gen+1)<SIGMA_STOP || F_opt_g(gen+1)<F_STOP);

    %% Save data
    num_gen(i) = gen+1;
    num_feval(i) = gen*LAM;
    r_g_t(:,i) = r_g;       % remove if memory demand gets too large; r_avg calcualtes avg. without saving all trials
    
    %% Check success (F_STOP, N_G)
    if F_opt_g(gen+1) < F_STOP || r_g(gen+1) < R_STOP
        r_avg = r_avg + r_g;  %r_g initialized with "nan", fastest run will set the mark by inserting nan; 
        f_avg = f_avg + F_opt_g;
        sigma_avg = sigma_avg + sigma_g;
        sigma_norm_avg = sigma_norm_avg + sigma_g./r_g*N;
        flag_success(i) = 1;
    else
        % flag_success: initialized with 0
    end
    
end % trials 
toc
eval = EvalDyn;
%% END: main loop

%% Fig
figure; hold on;
num_success = sum(flag_success);
title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES',', $\alpha$=$2\pi$', ', $A$=',num2str(A),', $N$=',num2str(N), ', ' sig_ctrl_str, ...
    ', $P_S$=',num2str(num_success/TRIALS, 2),' (',num2str(num_success),'/',num2str(TRIALS),')'];

%% Plot all dynamics
legend('AutoUpdate','off');
if PLOT_ALL_R==1
    for i=1:TRIALS
        if flag_success(i)==1, c = [0.6,0.6,0.6]; else, c = [0.9,0.9,0.9]; end
        plot(0:1:N_G, r_g_t(:,i), 'color', c);       % not normalized
    end
end
legend('AutoUpdate','on');

%% Averaging
if PLOT_AVG==1
    if sum(flag_success) == 0
        r_g_t_succ = r_g_t;
    else
        r_g_t_succ = r_g_t(:,logical(flag_success'));
    end
    g=0:1:N_G; 
        %plot(g, mean(r_g_t_succ, 2, 'includenan'), 'k-', 'DisplayName', 'mean($R$)');
        plot(g, median(r_g_t_succ, 2, 'includenan'), 'k-', 'DisplayName', 'median($R$)');
        %plot(g,  10.^mean(log10(r_g_t_succ), 2, 'includenan'), 'c-.', 'DisplayName', 'mean(log($R$))');
        %plot(g,  10.^median(log10(r_g_t_succ), 2, 'includenan'), 'm-.', 'DisplayName', 'test');
        %plot(g, mean(r_g_t_succ, 2, 'omitnan'), 'c--', 'DisplayName', 'mean($R(g)$), omitnan');
    dyn = EvalDyn;
%         dyn.dyn_avg_glob_attr(r_g_t, flag_success)
%         dyn.dyn_avg_btw_target(r_g_t, flag_success, 150, {'mean'})
%         dyn.dyn_avg_btw_target(r_g_t, flag_success, 150, {'median'})
end

%% FIG closing
title(title_str);
xlabel('$g$'), ylabel('$R(g)$');
set(gca, 'YScale', 'log');

fct_iter(MU, LAM,A, N, ALPHA, SIGMA_NORM, R_STOP, F_STOP, N_G, SEED, 'PHI1' )
fct_iter(MU, LAM,A, N, ALPHA, SIGMA_NORM, R_STOP, F_STOP, N_G, SEED, 'Y' )
myfigstyle(gcf, 8, 6, 9, 9);


%% Historgram
% dyn.histogram_gen(g, r_g_t, 0.8, 1.2)
% dyn.histogram_R(r_g_t, 75, 'lin')
