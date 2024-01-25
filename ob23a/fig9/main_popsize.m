
% clear

%% PARAMS
A_LIST = [0.1,1,4,7,10];
MU_RANGE = [1,10;  10,150;  30,350; 50,650; 100,850];
NUM_MU = 15;
TRIAL_LIST = [1000,1000,1000,1000,1000];
SEED = 2;

VARTHETA = 0.25;
N = 100;
N_G = 12000;
ALPHA = 2*pi;
TAU = 1/sqrt(2*N);
tau_str = '1/sqrt(2*N)';

SIGMA_NORM = nan; % nan for sigmaSA
SIGMA_STOP = 1e-8;
R_STOP = 1e-6;

%% R init
a = 100;   % R > 7*sqrt(N) for A=10 necessary
R_0 = a*sqrt(N);
rng(SEED); v = randn(N,1); 
Y_0 = v/norm(v)*R_0;
Y_HAT = zeros(N, 1);

%% Initialize variation
NUM_A = length(A_LIST);
% RES_MU_LIST = nan*zeros(NUM_MU,NUM_A);
% RES_PS = nan*zeros(NUM_MU, NUM_A);

tic
for k = NUM_A:-1:1

    %% SET POSITION
    A = A_LIST(k); 
    FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1);
    fprintf(' \t A = %i \n', A)
    
    %% GET MU(N);
    list_all = ceil(linspace(MU_RANGE(k,1), MU_RANGE(k,2), NUM_MU)');
    list_unq = unique(list_all);
    RES_MU_LIST(1:length(list_unq),k) = list_unq;
    
    TRIALS = TRIAL_LIST(k);
    
    for m = NUM_MU:-1:1
        
        MU = RES_MU_LIST(m,k); 
        if ~isnan(MU) && isnan(RES_PS(m,k))
            fprintf(' \t  MU = %i (%i/%i); (N: %i/%i) \n', MU, m, NUM_MU, k, NUM_A);
            LAM = ceil(MU/VARTHETA);

            %% Set MU-dep SIGMA_0
            SIGMA_0 = MU*e_vartheta_a_b(MU/LAM, 1, 0)*R_0/N;

            %% Initialize trial
            flag_success = zeros(TRIALS, 1);

            %% Main Loop 
            parfor i = 1:TRIALS

                if isnan(SIGMA_NORM)
                    [y_opt, F_opt_g, r_g, sigma_g, gen] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, N_G, 0);
                else
                    [y_opt, F_opt_g, r_g, sigma_g, gen] = muComLam_sigNorm(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_NORM, SIGMA_STOP, R_STOP, N_G, 0);
                end

                %% Increase max. value of N_G
                assert(gen+1<N_G)

                %% Check success (F_STOP, N_G)
                if r_g(gen+1) < R_STOP
                    flag_success(i) = 1;
                else
                    flag_success(i) = 0;
                end

            end % trials 

            RES_PS(m,k) = sum(flag_success)/TRIALS;
            save('wsp.mat');
        end % not isnan

    end % FOR: MU
end % FOR: N
toc

TRIAL_LIST = repmat(TRIAL_LIST,NUM_MU,1);

data = struct;
for n = 1:length(A_LIST)
    data(n).A = A_LIST(n);
    data(n).MU = RES_MU_LIST(:,n);
    data(n).TR = TRIAL_LIST(:,n);
    data(n).PS = RES_PS(:,n);
end
es.ALPHA = ALPHA;
es.N = N;
es.vartheta= VARTHETA;
es.tau = tau_str;
es.N_G = N_G;
es.a = a;
save('data_wsp.mat', 'data', 'es');