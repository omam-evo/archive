clear
addpath('../code');

%% CONFIG A
% A = 1;

%% CONFIG B
A = 7;

%% PARAMS
MU = 100;
LAM = 200; 
N = 100; 
TAU = 1/sqrt(8*N);

ALPHA = 2*pi;
G = 10000;
TRIALS = 100;
rng(2);

SIGMA_STOP = 1e-3;
R_STOP = 1e-2;
F_STOP = nan; %1e-2;
FEVAL_STOP = nan;

%% Plot options
PLOT_AVG = 1;
PLOT_ALL_R = 1;

%% Position fix
Y_0 = 50*ones(N, 1); 
R_0 = norm(Y_0);

%% Position random
% a = 20;   % R > 7*sqrt(N) for A=10 necessary
% R_0 = a*sqrt(N);
% SEED = 1; rng(SEED); v = randn(N,1); 
% Y_0 = v/norm(v)*R_0;

%% Set SIGMA
C_MU_LAM = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
C_VARTHETA = e_vartheta_a_b(MU/LAM, 1, 0);
x_2ndZero = sqrt( sqrt(N*(8*C_MU_LAM^2*MU^2 + N)) - N);
SIGMA_0 = x_2ndZero*R_0/N;

%% Fitness
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1); fit_name = 'Rastrigin';
Y_HAT = zeros(N, 1);

%% Get all extrema
ras_tools = RasTools;
[order_y_min, order_y_max] = ras_tools.get_all_extrema(ALPHA, A, N);

%% maxinfo
flag_success = zeros(TRIALS, 1);
r_g_t = zeros(G, TRIALS); 
sigma_g_t = zeros(G, TRIALS); 
% y_g_t = zeros(N, N_G+1, TRIALS); 

%% Main Loop 
tic
for i = 1:TRIALS

    fprintf('  t = %i \n', i)

    [y, f_g, r_g, sigma_g, gen, feval] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, FEVAL_STOP, G, 0);
   
    %% Save data
    num_gen(i) = gen;
    num_feval(i) = feval;
    r_g_t(:,i) = r_g;           % remove if memory demand gets too large; r_avg calcualtes avg. without saving all trials
    sigma_g_t(:,i) = sigma_g;
    %y_g_t(:,:,i) = y_g';       %  %maxinfo
    
    %% Check success (F_STOP, N_G)
    if f_g(gen) < F_STOP || r_g(gen) < R_STOP
        flag_success(i) = 1;
    else
        % flag_success: initialized with 0
    end
    
end % trials 
toc
%% END: main loop

P_S = sum(flag_success)/TRIALS;

%% GET MEDIANS
% select = flag_success;
% [~,~,~,fig1,fig2,fig3] = get_medians(r_g_t, sigma_g_t, flag_success, select, N, {'dyn','all','sign_g'});

%% A=1
if A==1
    select = flag_success;
    [~,~,~,fig1,fig2,fig3] = get_medians(r_g_t, sigma_g_t, flag_success, select, N, {'dyn','sign_r'});
    
    figure(fig1);
    ylim([1e-2,1e3]); 
    yticks([1e-2,1e-1,1,1e1,1e2,1e3]);
    xticks([0:200:1000]);
    myfigstyle(fig1, 6, 5, 8, 8);
    
    figure(fig2);
    ylim([20,60]); 
    yticks([20:10:60]);
    xticks([0:200:1000]);
    myfigstyle(fig2, 6, 5, 8, 8);
end

% A=7
if A==7
    select = ~flag_success;
    [~,~,~,fig1,fig2,fig3] = get_medians(r_g_t, sigma_g_t, flag_success, select, N, {'dyn','sign_r'});
    
    figure(fig1);
    ylim([1e-1,1e3]); 
    yticks([1e-1,1,1e1,1e2,1e3]);
    xticks([0:500:2000]);
    xlim([0,2000])
    myfigstyle(fig1, 6, 5, 8, 8);
    
    figure(fig2);
    ylim([0,60]); 
    yticks([0:20:60]);
    xticks([0:500:2000]);
    xlim([0,2000])
    myfigstyle(fig2, 6, 5, 8, 8);
end


% myfigstyle(fig2, 6, 5, 8, 8);

%% Fig
% title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES',', $\alpha$=$2\pi$', ', $A$=',num2str(A),', $N$=',num2str(N), ', $\tau$=' num2str(TAU), ...
%     ', $P_S$=',num2str(P_S,2)];

