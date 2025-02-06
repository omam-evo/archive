clear

%% Variation of N
% MODE_VAR = 1;
% N_LIST = [10,20,50,100,200,500,1000];
% MU_LIST = 1000*ones(1,length(N_LIST));

%% Variation of MU
MODE_VAR = 2;
MU_LIST = [10,30,100,300,1000,3000,10000];
N_LIST = 100*ones(1,length(MU_LIST));

%% Variation of CSA: (un-)comment sections (4a,4b,4c) below
% see below in section 'if strcmpi(MODE, 'CSA')'

%% General
TRIALS = 10;
G = 20000;
SHOW_ALL_DYN = 0;

RES = nan*zeros(length(N_LIST), TRIALS);
SIGN = nan*zeros(length(N_LIST),1);
COND_1 = nan*zeros(length(N_LIST),1);
COND_2 = nan*zeros(length(N_LIST),1);

%% Misc. init
MODE = 'CSA';
if strcmpi(MODE, 'CSA')
    S_0 = 0;
    %% CONFIG with CSA (4a)
    % MODE_CSA = 1;
    % E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    % C_handle = @(N,MU) 1/sqrt(N);
    % D_handle = @(N,MU) 1/C_handle(N,MU);
    %% CONFIG with CSA (4b)
    MODE_CSA = 1;
    E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    C_handle = @(N,MU) 1/N;
    D_handle = @(N,MU) 1/C_handle(N,MU);
    %% CONFIG with CSA (4c)
    % MODE_CSA = 2;
    % E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    % C_handle = @(N,MU) (MU+2)/(MU+N+5);
    % D_handle = @(N,MU) 1 + C_handle(N,MU) + 2 * max(0, sqrt((MU-1)/(N+1))-1);   
elseif strcmpi(MODE, 'SA')
    tau_handle = @(N) 1/sqrt(8*N); %1/sqrt(N)*sqrt((1-0.5^2)/(1-1/(2*0.8*sqrt(2*N)))) 
    SET_META_EP = 0;
end
sph = Sphere;

%% Position random
R_0 = 1e6;
rng(1); 

% figure; hold on;
for j=1:length(N_LIST)
    
    N = N_LIST(j);
    MU = MU_LIST(j);
    LAM = MU*2;

    FIT = Fitness('Sphere',N,[]);
    SIGNORM_ZERO = sqrt( sqrt(N*(8*e_vartheta_a_b(MU/LAM,1,0)^2*MU^2 + N)) - N);
    SIGMA_0 = SIGNORM_ZERO * R_0 / N;  

    SIGMA_STOP = 1e-6; %1e-10;
    R_STOP = 1;
    F_STOP = nan;

    %% Initialize  
    num_gen = zeros(TRIALS, 1);
    num_feval = zeros(TRIALS, 1);
    r_g_t = zeros(G, TRIALS); 
    sigma_g_t = zeros(G, TRIALS); 
    % y_g_t = zeros(N, G, TRIALS); 

    %% Main Loop 
    tic
    if SHOW_ALL_DYN==1
        figure; hold on;
    end
    for i = 1:TRIALS
        v = randn(N,1); 
        Y_0 = v/norm(v)*R_0;
        
        fprintf('  t = %i \n', i)

        if strcmpi(MODE, 'CSA')
            E_chi = E_chi_handle(N);
            C = C_handle(N,MU);
            D = D_handle(N,MU);
            [y, f_g, r_g, sigma_g, g_end, ~] = csa_es(FIT, N, MU, LAM, Y_0, SIGMA_0, S_0, C, D, E_chi, SIGMA_STOP, R_STOP, F_STOP, G, 0, MODE_CSA);
        elseif strcmpi(MODE, 'SA')
            TAU = tau_handle(N);
            [y, f_g, r_g, sigma_g, g_end, ~] = muComLam_sSA(SET_META_EP, FIT, N, MU, LAM, Y_0, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, G, 0);
        end
        if SHOW_ALL_DYN==1
            plot(r_g, 'k-');
            plot(sigma_g, 'r--');
        end

        %% Save data
        num_gen(i) = g_end;
        num_feval(i) = g_end*LAM;
        r_g_t(:,i) = r_g;           % remove if memory demand gets too large; r_avg calcualtes avg. without saving all trials
        sigma_g_t(:,i) = sigma_g;
        %y_g_t(:,:,i) = y_g';       %  %maxinfo

        RES(j,i) = g_end;
        if sigma_g(g_end)<SIGMA_STOP
            RES(j,i) = nan;
        end
        if TRIALS == 1
            plot(r_g)
            set(gca, 'YScale', 'log');
            legend;
        end
        
    end % trials 
    toc
    %% END: main loop
  
end

save('wsp.mat');

Gen = @(N,GAMMA) sqrt(N)./(1-GAMMA.^2)*log(R_0/R_STOP)/sqrt(2)/e_vartheta_a_b(MU/LAM,1,0);
if MODE_VAR ==1
    figure; hold on;
        x = linspace(min(N_LIST),max(N_LIST));
        base = 5;
        facs = base.^(0:1:2);
        % for p=1:length(facs)
        %     fac = facs(p);
        %     if p==1
        %         legend('autoupdate','on');
        %     else
        %         legend('autoupdate','off');
        %     end
        %     %plot(x, fac*x.^(1/2), 'c:', 'DisplayName', '$\sqrt{N}$');
        %     %plot(x, fac*x, 'm:', 'DisplayName', '$N$');
        % end
        if strcmpi(MODE, 'CSA')
            plot(x, Gen(x,0.8), 'k--');
            plot(x, Gen(x,0.9), 'k-.');
            plot(x, Gen(x,0.95), 'k--');
        elseif strcmpi(MODE, 'SA')
            plot(x, Gen(x,0.0), 'k--');
            plot(x, Gen(x,0.7), 'k--');
            plot(x, Gen(x,0.9), 'k--');
            plot(x, Gen(x,0.95), 'k--');
        end
        plot(N_LIST, mean(RES,2), 'b-');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('$N$');
        ylabel('Gen. $g$ with $R(g)/R(0)$=$10^{-6}$')
        myfigstyle(gcf,16,10,10,10);
else
    figure; hold on;
        %x = linspace(min(MU_LIST),max(MU_LIST));
        % base = 5;
        % facs = base.^(0:1:2);
        % for p=1:length(facs)
        %     fac = facs(p);
        %     if p==1
        %         legend('autoupdate','on');
        %     else
        %         legend('autoupdate','off');
        %     end
        %     plot(x, fac*x.^(1/2), 'c:', 'DisplayName', '$\sqrt{N}$');
        %     plot(x, fac*x, 'm:', 'DisplayName', '$N$');
        % end
        plot(MU_LIST, Gen(N_LIST,0.80), 'k--');
        plot(MU_LIST, Gen(N_LIST,0.9), 'k-.');
        plot(MU_LIST, Gen(N_LIST,0.95), 'k--');
        plot(MU_LIST, mean(RES,2), 'k-');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('$\mu$');
        ylabel('Gen. $g$ with $R(g)/R(0)$=$10^{-6}$')
        yscale('log');
        myfigstyle(gcf,16,10,10,10);
end