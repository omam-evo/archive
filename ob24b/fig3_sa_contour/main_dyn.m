clear all
addpath('../code');

%% Config A1
MODE = 'SA';
ALPHA = 2*pi;

%% Config A2
% MODE = 'SA';
% ALPHA = 10*pi;

%% Config B1
% MODE = 'METAEP';
% ALPHA = 2*pi;

%% Config B2
% MODE = 'METAEP';
% ALPHA = 10*pi;

%% General
N = 20;
A = 5;
FIT = Fitness('Rastrigin',N,[A, ALPHA]);
MU = 1000;
LAM = 2000;
TRIALS = 100;
verbose = 0;
G = 2000;

%% SET
if strcmpi(MODE, 'SA') || strcmpi(MODE, 'METAEP')
    TAU = 1/sqrt(2*N); 
    OPTIONS.weights = 'I';
elseif strcmpi(MODE, 'CSA')
    C = 1/sqrt(N);  %sqrt(N):fast, N:slow
    D = 1/C;
    E_chi = N; %sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));  %HGB: N
    S_0 = 0;
    OPTIONS.weights = 'I';
elseif strcmpi(MODE, 'CMA')
    OPTIONS.mode = 'on'; % 'on', 'off'
    OPTIONS.weights = 'W'; % 'W', 'I'
    OPTIONS.csa = 1; % 0:HGB, 1:HAN
    OPTIONS.S_0 = 0; 
    OPTIONS.P_0 = 0; 
    OPTIONS.cw = 1; % toggle rank-mu update 
    OPTIONS.matrix_norm = 0;
elseif strcmpi(MODE, 'SIGNORM')
    dummy = 1;
end

%% RNG
SEED = 1; 
rng(SEED);

%% Sphere
sph = Sphere;
SIGMA_NORM = sqrt( sqrt(N*(8*e_mu_lam_a_b_v2(MU,LAM,1,0)^2*MU^2 + N)) - N);
% [SIGMA_NORM, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});

%% Y via R
% R_0 = 100*sqrt(N);
% v = randn(N,1);
% Y_0 = v/norm(v)*R_0;
%% Y via R
Y_0 = 100*ones(N,1);

%% Sigma
% SIGMA_0 = SIGMA_NORM * norm(Y_0) / N;
SIGMA_0 = 1000;

%% Stop
SIGMA_STOP = 1e-10;
R_STOP = nan;
F_STOP = 1e-8; %1e-3;

%% Sphere
% FIT = Fitness('Sphere',N,[]);

%% Initialize
flag_success = nan*zeros(TRIALS, 1);
feval_t = zeros(TRIALS, 1);
r_g_t = zeros(G, TRIALS); 
sigma_g_t = zeros(G, TRIALS); 
y_g_t = zeros(N, G, TRIALS); 

%% Main Loop 
tic
for i = 1:TRIALS
    %v = randn(N,1); 
    %Y_0 = v/norm(v)*R_0;

    fprintf('  t = %i \n', i)
    if strcmpi(MODE, 'SA')
        [y, f_g, r_g, sigma_g, gen, y_g] = muComLam_sSA(FIT, N, MU, LAM, Y_0, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, G, verbose);
    elseif strcmpi(MODE, 'METAEP')
        [y, f_g, r_g, sigma_g, gen, y_g] = muComLam_metaEP(FIT, N, MU, LAM, Y_0, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, G, verbose);
    elseif strcmpi(MODE, 'CSA')
        [y, f_g, r_g, sigma_g, gen, y_g] = csa_es(FIT, N, MU, LAM, Y_0, SIGMA_0, S_0, C, D, E_chi, SIGMA_STOP, R_STOP, F_STOP, G, verbose);
    end

    %% Save data
    feval_t(i) = gen*LAM;
    r_g_t(:,i) = r_g;  
    sigma_g_t(:,i) = sigma_g;
    y_g_t(:,:,i) = y_g;
    
    if r_g(gen) < R_STOP || f_g(gen) < F_STOP, flag_success(i) = 1; else, flag_success(i) = 0; end
    
end % trials 
toc
%% END: main loop
P_S = sum(flag_success)/TRIALS;

%% MULTI EVAL: 'dyn','sign_g', 'all', 'sign_r'
select = ones(TRIALS,1); % select the single trial
[r_med,sigma_med,sigma_norm_med] = get_medians(r_g_t, sigma_g_t, flag_success, select, N, {'all', 'sign_r'}); %'sig_g','sign_g'$

%% Plot commands
% set(gca, 'XScale', 'log');
% myfigstyle(gcf,6,5,9,9)
% str = 'ln/ep_ras_colors'; savefig(gcf,[str, '.fig']); exportgraphics(gcf,[str, '.png'], 'resolution', 600, 'ContentType','vector');

%% Random colors
% L = findobj(gca,'Type','line');
% cc = rand(length(L),3);
% for i=1:length(L), L(i).Color = cc(i,:), end
