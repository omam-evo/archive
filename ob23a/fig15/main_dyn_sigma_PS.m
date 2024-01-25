clear
SAVE = 0;

%% CONFIG
MU = 400;
LAM = 800;
N = 100;
A = 10;
SIGMA_LIST = [0.1,0.2,0.3,(0.4:0.02:0.6),0.7,0.8,0.9,1];
N_G = 10000;
SEED = 1;

%% Additional config
ALPHA = 2*pi;
TRIALS = 500;

%% Position
a = 100;
R_0 = a*sqrt(N);
rng(SEED); 
v = randn(N,1); 
Y_0 = v/norm(v)*R_0;

%% CHECK if intersection btw. noise and sigma_0
% isolines_crit_point(MU, LAM, N, A, ALPHA, 101, 20, 1, 1, sigma_0)

%% Sphere and coefficients
%sph = Sphere;
%[x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});
%C_MU_LAM = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
%c_vartheta = e_vartheta_a_b(MU/LAM, 1, 0);

%% Fitness
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1);
Y_HAT = zeros(N, 1);

%% Set sigma(g) constant
sig_ctrl_str = ['$\sigma(g)$'];

%% SIGMA(G)
SIGMA_STOP = 1e-6;
R_STOP = 1e-4;
F_STOP = nan;

%% Main Loop 
flag_success = nan*zeros(length(SIGMA_LIST), TRIALS);
tic

for s=1:length(SIGMA_LIST)
    
    fprintf(' sigma=%f \n',SIGMA_LIST(s));
    
    if s==1 && TRIALS==1
        figure; hold on;
    end
    sigma_0 = SIGMA_LIST(s);
    
    %% Set constant
    SIGMA_G = sigma_0*ones(N_G+1,1);   
    
    %% Set sigma(g) decreasing
    gen = 1000;
    pow = log(sigma_0/SIGMA_STOP)/log(10);
    sigma_ramp = sigma_0*10.^linspace(0, -(pow+1e-3), gen+1)';
    SIGMA_G(end-gen:end) = sigma_ramp;

    %% Initialize
    init_vec = zeros(N_G+1, 1); % first value contains 0-th generation
    init_trial = zeros(TRIALS, 1);
    [sigma_avg, sigma_norm_avg] = deal(init_vec, init_vec);
    %[r_avg, f_avg] = deal(init_vec, init_vec);
    %[num_gen, num_feval] = deal(init_trial, init_trial);

    %% maxinfo
    %r_g_t = zeros(N_G+1, TRIALS); 
    %sigma_g_t = zeros(N_G+1, TRIALS); 
    %y_g_t = zeros(N, N_G+1, TRIALS); 

    for i = 1:TRIALS
            fprintf(' \t trial=%i/%i \n',i,TRIALS);

            [y_opt, F_opt_g, r_g, sigma_g, gen, y_g] = muComLam_sigma_maxInfo(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_G, SIGMA_STOP, R_STOP, N_G, 0);

            %% Save data
            %num_gen(i) = gen+1;
            %num_feval(i) = gen*LAM;
            %r_g_t(:,i) = r_g;           % remove if memory demand gets too large; r_avg calcualtes avg. without saving all trials
            %sigma_g_t(:,i) = sigma_g;
            %y_g_t(:,:,i) = y_g';       %  %maxinfo

            %% Check success (F_STOP, N_G)
            if F_opt_g(gen+1) < F_STOP || r_g(gen+1) < R_STOP
                %r_avg = r_avg + r_g;  %r_g initialized with "nan", fastest run will set the mark by inserting nan; 
                %f_avg = f_avg + F_opt_g;
                %sigma_avg = sigma_avg + sigma_g;
                %sigma_norm_avg = sigma_norm_avg + sigma_g./r_g*N;
                flag_success(s,i) = 1;
            else
                flag_success(s,i) = 0;
            end

    end % trials 
    save('wsp.mat')
    if TRIALS==1
        R_g = r_g_t;
        gen = 0:N_G; gen(isnan(R_g))=nan; 
        xlabel('$g$'); ylabel('$R(g)$-dynamics');
        plot(gen, R_g, '-', 'DisplayName', ['$\sigma=$',num2str(sigma_0)]);
        set(gca, 'YScale', 'log');
    end
end % sigma
toc

if TRIALS==1
    legend;
    title_str = ['Rastrigin: (',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha$=$2\pi$', ', $A$=',num2str(A), ', $N$=',num2str(N)];
    title(title_str);
    myfigstyle(gcf, 16, 10, 10, 10);

    gen = 1000;
    figure;
        histogram(y_g_t(:,end-(gen-1):end), 'BinWidth', 0.02);
        xlabel('$y_i$ ($i=1,...,N$)');
        ylabel(['Counts of $y_i$ (last ', num2str(gen) ' gen)' ]);
        title(title_str);
    myfigstyle(gcf, 16, 10, 10, 10);
else
    P_S = sum(flag_success,2)/TRIALS;
    figure; hold on;
    plot(SIGMA_LIST,P_S, '.')
    xlabel('$\sigma$');
    ylabel('$P_S$');
    %title_str = ['Rastrigin: (',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha$=$2\pi$', ', $A$=',num2str(A), ', $N$=',num2str(N)];
    %title(title_str);
    xlim([SIGMA_LIST(1),SIGMA_LIST(end)]);
    xticks([SIGMA_LIST(1):0.1:SIGMA_LIST(end)]);
    yticks(0:0.2:1);
    myfigstyle(gcf, 8, 4.5, 9, 9); 
end