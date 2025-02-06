clear

% MU_LIST = [10,30,100,300,1000,3000,10000];
MU_LIST = [10,100,1000];

%% CONFIG ABC1
N_LIST = 50*ones(1,length(MU_LIST));
TRIALS = 1;
colors = {'r','b','k'};

VARTHETA = 0.5;

verbose = 0;
G = 20000;
g0 = 50;

%% Misc
MODE = 'CSA';
if strcmpi(MODE, 'CSA')
    S_0 = 0;
    %% AR, HGB derandomized 1/sqrt
    % MODE_CSA = 0;
    % E_chi_handle = @(N) N;
    % C_handle = @ (N,MU) 1/sqrt(N);
    % D_handle = @(N,MU) 1/C_handle(N,MU);
    %% CONFIG A: Hansen derandomized 1/sqrt
    % MODE_CSA = 1;
    % E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    % C_handle = @(N,MU) 1/sqrt(N);
    % D_handle = @(N,MU) 1/C_handle(N,MU);
    %% CONFIG A2: Hansen derandomized 1/sqrt
    % MODE_CSA = 1;
    % E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    % C_handle = @(N,MU) 1/N;
    % D_handle = @(N,MU) 1/C_handle(N,MU);
    %% CONFIG B: Hansen tutorial default
    MODE_CSA = 2;
    E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    C_handle = @(N,MU) (MU+2)/(MU+N+5);
    D_handle = @(N,MU) 1 + C_handle(N,MU) + 2 * max(0, sqrt((MU-1)/(N+1))-1);
    %% Hansen tutorial modified (mimic derandomized 1/sqrt)
    % MODE_CSA = 2;
    % E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    % C_handle = @(N,MU) 1/sqrt(N);
    % D_handle = @(N,MU) 1 + C_handle(N,MU);   
elseif strcmpi(MODE, 'SA')
    tau_handle = @(N) 1/sqrt(2*N); %1/sqrt(N)*sqrt((1-0.5^2)/(1-1/(2*0.8*sqrt(2*N)))) 
end


%% Y via R
R_0 = 10^3;
SIGMA_STOP = 1e-15;
R_STOP = 1e-6;
F_STOP = nan; %1e-3;
% SIGMA_0: initilaized as fct of MU

%% Initialize

flag_success = nan*zeros(TRIALS, 1);
feval_t = zeros(TRIALS, 1);
r_g_t = zeros(G, TRIALS); 
f_g_t = zeros(G, TRIALS); 
sigma_g_t = zeros(G, TRIALS); 
signorm_g_t = zeros(G, TRIALS); 
%y_g_t = zeros(N, G, TRIALS); 
sign_ss_per_trial = nan*zeros(length(MU_LIST),TRIALS);

sph = Sphere;
mode = 'full';
sig_ss_res = nan*zeros(length(MU_LIST),1);
signzero_approx = nan*sig_ss_res;
sign_phi0_full = nan*sig_ss_res;
sign_opt = nan*sig_ss_res;
generations = zeros(length(MU_LIST),TRIALS);

%% Main Loop 
tic
fig_dyn = figure;

for m=1:length(MU_LIST)
    
    MU = MU_LIST(m);
    N = N_LIST(m);
    LAM = round(MU/VARTHETA);
    fprintf('mu = %i, N = %i \n', MU, N)
    
    FIT = Fitness('Sphere',N,[]);
    Y_0 = 1*ones(N,1);
    SIGNORM_ZERO = (8*N)^(1/4)*sqrt(e_vartheta_a_b(MU/LAM,1,0)*MU);
    SIGMA_0 = 10;  

    rng(1);
    for i=1:TRIALS
        fprintf('  t = %i \n', i)
        if strcmpi(MODE, 'SA')
            TAU = tau_handle(N);
            [y, f_g, r_g, sigma_g, g_end, ~] = muComLam_sSA(FIT, N, MU, LAM, Y_0, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, G, verbose);
        elseif strcmpi(MODE, 'CSA')
            E_chi = E_chi_handle(N);
            C = C_handle(N,MU);
            D = D_handle(N,MU);
            [y, f_g, r_g, sigma_g, g_end, ~] = csa_es(FIT, N, MU, LAM, Y_0, SIGMA_0, S_0, C, D, E_chi, SIGMA_STOP, R_STOP, F_STOP, G, verbose, MODE_CSA);
        end

        f_g_t(:,i) = f_g;
        sigma_g_t(:,i) = sigma_g;
        r_g_t(:,i) = r_g;
        signorm = sigma_g./r_g*N;
        signorm_g_t(:,i) = signorm;
        generations(m,i) = g_end;
        % sign_ss_per_trial(m,i) = mean(signorm(100:end), 'omitnan'); 

    end

    %% Evaluation from experiment
    figure(fig_dyn); hold on;

    plot(1:G, r_g_t, '-', 'color',colors{m},'DisplayName',['mu=',num2str(MU)]);
    plot(1:G, sigma_g_t, '-.', 'color',colors{m},'DisplayName',['mu=',num2str(MU)]);
    plot(1:G, signorm_g_t, ':', 'color',colors{m},'DisplayName',['mu=',num2str(MU)]);

    %signorm = sign_g_med(~isnan(sign_g_med)); % select all not-nan values 
    % sig_ss_res(m) = median(sign_g_med(g0:end),'omitnan');
    signorm = signorm_g_t(g0:end,:);
    signorm = signorm(~isnan(signorm));
    sig_ss_res(m) = median(signorm);

    %% N-dep progress rate formula
    [x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU,LAM,1,0), MU, LAM, N, {mode});
    sign_opt(m) = x_opt;
    sign_phi0_full(m) = x_2ndZero;

    %% Approximation
    signzero_approx(m) = (8*N)^(1/4)*(e_vartheta_a_b(MU/LAM,1,0)*MU)^(1/2);

    %% Additional stuff
    figure(fig_dyn); hold on;
    xlabel('$g$');
    ylabel('Dynamics')
    set(gca, 'YScale', 'log');
    myfigsize(fig_dyn,8,4,9,9);
    % yticks([1e-6,1e-3,1,1e3,1e6]);
    % ylim([1e-3,1e7]);
    
end % MU
toc

