clear

%% CONFIG 1
% N = 10;

%% CONFIG 2
N = 100;

%% Settings
MU = 5000;
LAM = 10000;
TRIALS = 1;
verbose = 0;
G = 2;
MODE = 'SA';
TAU = 0;

%% RNG
SEED = 1; 
rng(SEED);
% randn("state", 1);

%% Y via R
% SIGMA_NORM = 0.5;
R_0 = 10;
v = randn(N,1);
Y_0 = v/norm(v)*R_0;
%% Y via R
% Y_0 = 100*ones(N,1);

%% Sigma
SIGMA_0 = 1; %1;

%% Stop
SIGMA_STOP = 1e-2;
R_STOP = nan;
F_STOP = 1e-3; %1e-3;

%% Sphere
A = 10;
ALPHA = 2*pi;
FIT = Fitness('Rastrigin',N,[A,ALPHA]);

%% Initialize
flag_success = nan*zeros(TRIALS, 1);
feval_t = zeros(TRIALS, 1);
r_g_t = zeros(G, TRIALS); 
sigma_g_t = zeros(G, TRIALS); 
y_g_t = zeros(N, G, TRIALS); 

%% Main Loop 
tic
[y, f_g, r_g, sigma_g, gen, y_g] = muComLam_sSA(FIT, A, ALPHA, N, MU, LAM, Y_0, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, G, verbose);



