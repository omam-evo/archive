clear


%% CONFIG 1
% tau_str = '1/sqrt(N)'; 

%% CONFIG 2
% tau_str = '1/sqrt(2*N)'; 

%% CONFIG 3
tau_str = '1/sqrt(8*N)'; 

MU = 400;
LAM = 800;
N = 100;
A = 10;

ALPHA = 2*pi;
TRIALS = 100;
N_G = 5000;

sph = Sphere;
[x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});
C_MU_LAM = e_mu_lam_a_b_v2(MU, LAM, 1, 0);

%% Position random
R_0 = 100*sqrt(N);
SEED = 1; 
rng(SEED); 
v = randn(N,1); 
Y_0 = v/norm(v)*R_0;

%% SIGMA(G)
% MODE = 'SIGMA';
% %SIGMA_G = linspace(x_opt*R_0/N, 0, N_G+1); %linear
% SIGMA_G = C_MU_LAM*MU*R_0/N*2.^linspace(0, -20, N_G+1); 

%% SIGMA NORM
% MODE = 'SIGNORM';
% SIGMA_NORM = nan;

%% sigmaSA
MODE = 'SA';

SIGMA_STOP = 1e-5; %1e-10;
R_STOP = 1e-3;
F_STOP = nan;

%% Set SIGMA & TAU
if strcmpi(MODE, 'SA')
    SIGMA_0 = x_opt*R_0/N;
    sig_ctrl_str = ['$\tau$=', tau_str, ''];
    TAU = eval(tau_str);
elseif strcmpi(MODE, 'SIGNORM')
    sig_ctrl_str = ['$\sigma^*$=', num2str(SIGMA_NORM)];
elseif strcmpi(MODE, 'SIGMA')
    sig_ctrl_str = ['$\sigma(g)$'];
end

%% Fitness
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1);
Y_HAT = zeros(N, 1);

%% Initialize
init_vec = zeros(N_G+1, 1); % first value contains 0-th generation
init_trial = zeros(TRIALS, 1);
[flag_success, num_gen, num_feval] = deal(init_trial, init_trial, init_trial);

%% maxinfo
r_g_t = zeros(N_G+1, TRIALS); 
sigma_g_t = zeros(N_G+1, TRIALS); 
% y_g_t = zeros(N, N_G+1, TRIALS); 

%% Main Loop 
tic
for i = 1:TRIALS

    fprintf('  t = %i \n', i)
    [y_opt, F_opt_g, r_g, sigma_g, gen] = deal(0,0,0,0,0); %dummy init
    
    if strcmpi(MODE, 'SA')
        [y_opt, F_opt_g, r_g, sigma_g, gen] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, N_G, 0);
    elseif strcmpi(MODE, 'SIGNORM')
        [y_opt, F_opt_g, r_g, sigma_g, gen] = muComLam_sigNorm(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_NORM, SIGMA_STOP, R_STOP, N_G, 0);
	elseif strcmpi(MODE, 'SIGMA')
        [y_opt, F_opt_g, r_g, sigma_g, gen, ~] = muComLam_sigma_maxInfo(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_G, nan, R_STOP, N_G, 0);
    end

    %% Save data
    num_gen(i) = gen+1;
    num_feval(i) = gen*LAM;
    r_g_t(:,i) = r_g;           % remove if memory demand gets too large; r_avg calcualtes avg. without saving all trials
    sigma_g_t(:,i) = sigma_g;
    %y_g_t(:,:,i) = y_g';       %  %maxinfo
    
    %% Check success (F_STOP, N_G)
    if F_opt_g(gen+1) < F_STOP || r_g(gen+1) < R_STOP
        flag_success(i) = 1;
    else
        % flag_success: initialized with 0
    end
    
end % trials 
toc
%% END: main loop

% save('wsp.mat');
% 
% select = ~flag_success;
% P_S = sum(flag_success)/TRIALS
% get_medians(r_g_t, sigma_g_t, flag_success, select, N);
% 
% myfigstyle(gcf, 8, 5, 9, 9);
% xticks([0:20:120]);
% yticks([0:2:10]);
