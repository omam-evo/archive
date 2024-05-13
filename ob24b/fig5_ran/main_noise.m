clear
addpath('../code');

%% CONFIG A: NOISE
% SET_metaEP = 0;

%% CONFIG B: NOISE
SET_metaEP = 1;

%% PARAMS
N = 100;
FIT = @(x) 1*randn(1,size(x,2));
SIGMA_0 = 1;
DIST = 0;
MU = 100;
LAM = 200; 
TAU = 1/sqrt(2*N);

TRIALS = 100;
G = 10000;
rng(2);

SIGMA_STOP = 1e-6;
R_STOP = nan;
F_STOP = nan; %1e-2;
FEVAL_STOP = nan;

%% Plot options
PLOT_AVG = 1;
PLOT_ALL_R = 1;

%% Position fix
Y_0 = ones(N,1); 
% R_0 = norm(Y_0);

%% Position random
% a = 20;   % R > 7*sqrt(N) for A=10 necessary
% R_0 = a*sqrt(N);
% SEED = 1; rng(SEED); v = randn(N,1); 
% Y_0 = v/norm(v)*R_0;

%% Set SIGMA
C_MU_LAM = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
C_VARTHETA = e_vartheta_a_b(MU/LAM, 1, 0);
x_2ndZero = sqrt( sqrt(N*(8*C_MU_LAM^2*MU^2 + N)) - N);
Y_HAT = zeros(N, 1);

%% maxinfo
flag_success = zeros(TRIALS, 1);
r_g_t = zeros(G, TRIALS); 
sigma_g_t = zeros(G, TRIALS); 
y_g_t = zeros(N, G, TRIALS); 

%% Main Loop 
tic
for i = 1:TRIALS

    fprintf('  t = %i \n', i)
    if DIST==0
        [y, f_g, r_g, sigma_g, gen, feval, y_g] = muComLam_metaEP_noFclass(SET_metaEP, FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, FEVAL_STOP, G, 0);
    else
        [y, f_g, r_g, sigma_g, gen, feval, y_g] = muComLam_metaEP_dist_noFclass(SET_metaEP, FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, FEVAL_STOP, G, 0);
    end
    %% Save data
    num_gen(i) = gen;
    num_feval(i) = feval;
    r_g_t(:,i) = r_g;           % remove if memory demand gets too large; r_avg calcualtes avg. without saving all trials
    sigma_g_t(:,i) = sigma_g;
    y_g_t(:,:,i) = y_g;       %  %maxinfo
    
    %% Check success (F_STOP, N_G)
    if f_g(gen) < F_STOP || r_g(gen) < R_STOP
        flag_success(i) = 1;
    else
        % flag_success: initialized with 0
    end
    
end % trials 
toc
%% END: main loop

%% GET MEDIANS
select = flag_success;
get_means(r_g_t, sigma_g_t, flag_success, ~select, N, {'all','dyn','sig_g'});
% get_medians(r_g_t, sigma_g_t, flag_success, ~select, N, {'all','dyn','sig_g'});


%% sigma-growth lognormal
if SET_metaEP==0
    x=1:G;
    hold on; plot(x,exp(TAU^2/2*x),'c--');
    set(gca, 'YScale', 'log');
else
    hold on; plot([1,G],[1,1],'c--')
end
% myfigstyle(gcf, 30, 10, 10, 10);

%% DISTRIBUTION
% s = SIGMA_0*TAU;
% x = SIGMA_0 + linspace(-5*s,5*s,1001);
% figure; plot(x, 1/s*normpdf((x-SIGMA_0)/s));

%% AVERAGING PLOT
% figure; hold on;
% num_success = sum(flag_success);

%% Single Dynamics
% plot(r_g, 'k-');
% plot(sigma_g, 'r-');
% plot(sigma_g./r_g*N, 'b-');

%% Component Dynamics
% for i=1:N
%     plot(y_g(i,:).^2)
% end


%% HISTOGRAM
% figure; hold on;
%     title(title_str);
%     histogram(sigma_g_t(10000,:));


%% Plot all dynamics
% figure; hold on;
% if PLOT_ALL_R==1
%     for i=1:TRIALS
%         if flag_success(i)==1, c = [0.6,0.6,0.6]; else, c = [0.9290, 0.6940, 0.1250]; end
%         plot(0:1:G-1, r_g_t(:,i), 'color', c);       % not normalized
%     end
% end
% 
% r_g_t_succ = r_g_t(:,logical(flag_success'));
% r_g_t_unsucc = r_g_t(:,logical(~flag_success'));
% 
% plot(0:1:G-1, median(r_g_t_succ, 2, 'includenan'), 'k-');
% plot(0:1:G-1, median(r_g_t_unsucc, 2, 'includenan'), 'r-');
% 
% title(title_str);
% xlabel('$g$'), ylabel('$R(g)$');
% set(gca, 'YScale', 'log');
% myfigstyle(gcf, 16, 10, 10, 10);
