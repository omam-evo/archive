clear
addpath('../n0_generic');

%% CONFIG 1
MU = 100;
LAM = 200; 
N = 100; 
A = 1;
ALPHA = 2*pi;

%% CONFIG 2
% MU = 500;
% LAM = 2000; 
% N = 50; 
% A = 10;
% ALPHA = 2*pi;

%% Variation of TAU
% TAU = 1/sqrt(8*N);
% TAU = sqrt(1/N.*(1-8./(ALPHA.^2.*A).*lambertw(0, ALPHA.^3.*A.^(3/2)/8)));
TAU = 1/sqrt(N);

G = 5000;
TRIALS = 100;
rng(2);

SIGMA_STOP = 1e-5;
R_STOP = 1e-4;
F_STOP = nan; %1e-2;
FEVAL_STOP = nan;

%% Plot options
PLOT_AVG = 1;
PLOT_ALL_R = 1;

%% Position fix
Y_0 = 100*ones(N, 1); 
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
GEN_MEAN = mean(num_gen(logical(flag_success)));

%% GET MEDIANS
select = flag_success;

% SHOW MEDIAN SUCCESSFUL
get_medians(r_g_t, sigma_g_t, flag_success, select, N, {'sign_r'});

% SHOW MEDIAN UNSUCCESSFUL
get_medians(r_g_t, sigma_g_t, flag_success, ~select, N, {'sign_r'});

%% Fig
% figure; hold on;
% num_success = sum(flag_success);
% title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES',', $\alpha$=$2\pi$', ', $A$=',num2str(A),', $N$=',num2str(N), ', $\tau$=' num2str(TAU), ...
%     ', $P_S$=',num2str(P_S,2)];
% 
% %% Plot all dynamics
% legend('AutoUpdate','off');
% if PLOT_ALL_R==1
%     for i=1:TRIALS
%         if flag_success(i)==1, c = [0.6,0.6,0.6]; else, c = [0.9290, 0.6940, 0.1250]; end
%         plot(0:1:G-1, r_g_t(:,i), 'color', c);       % not normalized
%     end
% end
% legend('AutoUpdate','on');
% 
% %% Averaging
% if PLOT_AVG==1
% 
%     r_g_t_succ = r_g_t(:,logical(flag_success'));
%     r_g_t_unsucc = r_g_t(:,logical(~flag_success'));
% 
%     plot(0:1:G-1, median(r_g_t_succ, 2, 'includenan'), 'k-', 'DisplayName', 'median($R$)');
%     plot(0:1:G-1, median(r_g_t_unsucc, 2, 'includenan'), 'r-', 'DisplayName', 'median($R$)');
% end
% 
% %% FIG closing
% title(title_str);
% xlabel('$g$'), ylabel('$R(g)$');
% set(gca, 'YScale', 'log');
% myfigstyle(gcf, 16, 10, 10, 10);

