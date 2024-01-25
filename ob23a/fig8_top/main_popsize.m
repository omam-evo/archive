clear

N_LIST = [30,50,100,300,1000,3000];
MU_RANGE = [40,60; 55,100; 85,150; 180,250 ; 420,500; 700,1000];
RANGE_MODE = 'LIN';
NUM_MU = 5;
TRIAL_LIST = repmat([2000,2000,1000,700,500,400], NUM_MU,1); %[2000,2000,1000,700,500,400]
rng('shuffle');

VARTHETA = 0.25;
A = 1;
N_G = 5000;
ALPHA = 2*pi;
get_tau = @(N) 1/sqrt(2*N);
tau_str = '1/sqrt(2*N)';

SIGMA_NORM = nan; % nan for sigmaSA
SIGMA_STOP = 1e-5;
R_STOP = 1e-3;

%% R init
a = 10; % A=1 

%% Fitness
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1); fit_name = 'Rastrigin';

%% Initialize variation
NUM_N = length(N_LIST);
RES_MU_LIST = nan*zeros(NUM_MU,NUM_N);
RES_PS = nan*zeros(NUM_MU, NUM_N);

%% GET MU
for n = 1:NUM_N
    if strcmpi(RANGE_MODE, 'LOG')
        logfac = log(MU_RANGE(n,2)/MU_RANGE(n,1))/NUM_MU;
        fac = exp(logfac);
        list_all = round(MU_RANGE(n,1)*fac.^linspace(0, NUM_MU, NUM_MU)');
        %list_unq = unique(list_all);
        RES_MU_LIST(:,n) = list_all;
    elseif strcmpi(RANGE_MODE, 'LIN')
        list_all = round(linspace(MU_RANGE(n,1), MU_RANGE(n,2), NUM_MU)');
        %list_unq = unique(list_all);
        RES_MU_LIST(:,n) = list_all;
    elseif strcmpi(RANGE_MODE, 'MAN')
        RES_MU_LIST(:,n) = MU_RANGE(n,:);
    end
end

tic
for n = NUM_N:-1:1

    %% SET POSITION
    N = N_LIST(n); 
    fprintf(' \t N = %i \n', N)
    R_0 = a*sqrt(N);
    v = randn(N,1); 
    Y_0 = v/norm(v)*R_0;
    Y_HAT = zeros(N, 1); 
    
    %% Set N-dep TAU
    if isnan(SIGMA_NORM)
        TAU = get_tau(N);
    end
    
    for m = NUM_MU:-1:1
        
        TRIALS = TRIAL_LIST(m,n);
        MU = RES_MU_LIST(m,n); 
        if isnan(RES_PS(m,n))
            fprintf(' \t  MU = %i (%i/%i); (N: %i/%i) \n', MU, m, NUM_MU, n, NUM_N);
            LAM = ceil(MU/VARTHETA);

            %% Set MU-dep SIGMA_0
            sph = Sphere;
            [x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});
            SIGMA_0 = x_opt*R_0/N;

            %% Initialize trial
            flag_success = nan*zeros(TRIALS, 1);

            %% Main Loop 
            for i = 1:TRIALS
                if isnan(flag_success(i))
                    if isnan(SIGMA_NORM)
                        [y_opt, F_opt_g, r_g, sigma_g, gen] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, N_G, 0);
                    else
                        [y_opt, F_opt_g, r_g, sigma_g, gen] = muComLam_sigNorm(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_NORM, SIGMA_STOP, R_STOP, N_G, 0);
                    end

                    %% Increase max. value of N_G
                    assert(gen+1<N_G)

                    %% Check success (F_STOP, N_G)
                    if F_opt_g(gen+1) < R_STOP
                        flag_success(i) = 1;
                    else
                        flag_success(i) = 0;
                    end
                end

            end % trials 

            RES_PS(m,n) = sum(flag_success)/TRIALS;
            save('wsp.mat');
        end % not isnan

    end % FOR: MU
end % FOR: N
toc

data = struct;
for n = 1:length(N_LIST)
    data(n).N = N_LIST(n);
    data(n).MU = RES_MU_LIST(:,n);
    data(n).TR = TRIAL_LIST(:,n);
    data(n).PS = RES_PS(:,n);
end
es.A = A;
es.ALPHA = ALPHA;
es.vartheta= VARTHETA;
es.tau = tau_str;
es.N_G = N_G;
es.a = a;
save('data_wsp.mat', 'data', 'es');

    