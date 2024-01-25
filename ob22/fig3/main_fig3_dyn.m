clear

%% CONFIG 1
MU = 150;
LAM = 300;

%% CONFIG 2
% MU = 1500;
% LAM = 3000;

%% PARAMS
A = 1;
N = 100;
TRIALS = 100;
N_G = 2000;
ALPHA = 2*pi;
SIGMA_0 = 1;
SIGMA_NORM = 30;
SIGMA_STOP = 1e-6;
F_STOP = 1e-6;

%% Position
Y_0 = 10;
% R_0 = sqrt(N*10^2);
% rng(1); v = randn(N,1); Y_0 = v/norm(v)*R_0;

%% Set TAU
% if MU == 1, TAU = 1/sqrt(N_DIM); else, TAU = 1/sqrt(2*N_DIM); end
if isnan(SIGMA_NORM)
    TAU = 1/sqrt(8*N);
    sig_ctrl_str = ['$\tau$=', num2str(TAU,2)];
else
    sig_ctrl_str = ['$\sigma^*$=', num2str(SIGMA_NORM)];
end

%% Fitness
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1); fit_name = 'Rastrigin';
Y_HAT = zeros(N, 1);

%% Get all extrema
% ras_tools = RasTools;
% [order_y_min, order_y_max] = ras_tools.get_all_extrema(ALPHA, A, N);

%% Settings
title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $A$=',num2str(A),', $N$=',num2str(N), ', ' sig_ctrl_str];

%% Initialize
init_vec = zeros(N_G+1, 1); % first value contains 0-th generation
init_trial = zeros(TRIALS, 1);
[sigma_avg, sigma_norm_avg] = deal(init_vec, init_vec);
[r_avg, f_avg] = deal(init_vec, init_vec);
[flag_success, num_gen, num_feval] = deal(init_trial, init_trial, init_trial);

%% Main Loop 
tic
parfor i = 1:TRIALS

    fprintf('  t = %i \n', i)
    rstream = RandStream('threefry4x64_20');
    rstream.Substream = i;
    
    if isnan(SIGMA_NORM)
        [y_opt, F_opt_g, delta_r_g, sigma_g, gen] = muComLam_sigSA_randStream(rstream, FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, F_STOP, N_G, 0);
    else
        [y_opt, F_opt_g, delta_r_g, sigma_g, gen] = muComLam_sigNorm_randStream(rstream, FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_NORM, SIGMA_STOP, F_STOP, N_G, 0);
    end
    
    %% Assertion failes if N_G is reached before F_STOP or SIGMA_STOP => increase max. value of N_G
    %assert(sigma_g(gen+1)<SIGMA_STOP || F_opt_g(gen+1)<F_STOP);

    %% Save data
    num_gen(i) = gen+1;
    num_feval(i) = gen*LAM;
    
    %% Check success (F_STOP, N_G)
    if F_opt_g(gen+1) < F_STOP && (gen+1)<N_G
        r_avg = r_avg + delta_r_g;
        f_avg = f_avg + F_opt_g;
        sigma_avg = sigma_avg + sigma_g;
        sigma_norm_avg = sigma_norm_avg + sigma_g./delta_r_g*N;
        flag_success(i) = 1;
    else
        % flag_success: initialized with 0
    end

end % trials 
toc

%% EVALUATE
num_success = sum(flag_success); % count successful runs
r_avg(r_avg==0) = nan; % r_avg initialized with 0 (values added during sim.); "0" no values were written => "nan"
mask_nan = ones(N_G+1,1); mask_nan(isnan(r_avg)) = nan; % mask_nan sets not written elements as "nan" (N_G >(!) largest occuring generation)
title_str = [title_str, ', $P_S$=',num2str(num_success/TRIALS, 2)]; % evaluate P_s
fprintf('EVAL: min(num_gen)=%i, max(num_gen)=%i \n', min(num_gen), max(num_gen));

% Average values over num_success
r_avg = r_avg/num_success .* mask_nan;
f_avg = f_avg/num_success .* mask_nan;
sigma_avg = sigma_avg/num_success .* mask_nan;
sigma_norm_avg = sigma_norm_avg/num_success .* mask_nan;

fig = figure;
    g=0:1:N_G; 
    hold on; title(title_str);
    xlabel('$g$');
    ylabel('Dynamics (avg.)');
    set(gca, 'DefaultLineLineWidth', 2);
        plot(g, r_avg, 'k-', 'DisplayName', '$R$');
        %plot(g, f_avg, 'DisplayName', '$f$');
        %plot(g, sigma_avg, 'DisplayName', '$\sigma$');
        %plot(g, sigma_norm_avg, 'DisplayName', '$\sigma^*$');
    set(gca, 'YScale', 'log');
    legend('location', 'best');
    yticks([1e-6, 1e-4, 1e-2, 1, 1e2, 1e4]);
    myfigstyle(fig, 12, 6, 9, 9);