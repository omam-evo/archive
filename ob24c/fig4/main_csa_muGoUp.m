clear

%% SET SIGMA RESCALE
LIST_SIGMA_RESCALE = [0, 0.5, 1]; %0: (25a); 0.5: (25b); 1: (25c)
%% SET Dimension
N = 10;  %{10,1000}
%% SET wait (0: no; 1: yes)
SET_WAIT = 0; %{0,1}
%% SET CSA variant below (4a,4b,4c)
% see 'if strcmpi(MODE, 'CSA')'

alpha_mu = 2;
G = 2000;
MU_MIN = 4;
MU_MAX = 1024;
THETA = 0.5;
TRIALS = 1;
VARTHETA = 0.5;
verbose = 0;

%% SWEEP MU
if SET_WAIT==0
    MU_SWEEP = MU_MIN * ones(G,1);
    G_0 = 200;
    up_down = 1;
    mu = MU_MIN;
    for i=G_0:G
        MU_SWEEP(i) = mu;
        if up_down == 1
            mu = ceil(mu*alpha_mu);
            if mu>MU_MAX
                mu = MU_MAX;
                up_down = -1;
            end
        else
            mu = floor(mu/alpha_mu);
            if mu<MU_MIN
                mu = MU_MIN;
                up_down = 1;
            end
        end
    end
    MU_MED = median(MU_SWEEP(G_0:end));
end

%% SWEEP MU with WAIT
if SET_WAIT==1
    L = ceil(sqrt(N));
    MU_SWEEP = MU_MIN * ones(G,1);
    G_0 = 200;
    up_down = 1;
    mu = MU_MIN;
    for i=G_0:L:G
        MU_SWEEP(i:i+L) = mu;
        if up_down == 1
            mu = ceil(mu*alpha_mu);
            if mu>=MU_MAX
                mu = MU_MAX;
                up_down = -1;
            end
        else
            mu = floor(mu/alpha_mu);
            if mu<=MU_MIN
                mu = MU_MIN;
                up_down = 1;
            end
        end
    end
    MU_MED = median(MU_SWEEP(G_0:end));
end



%% Misc
MODE = 'CSA';
if strcmpi(MODE, 'CSA')
    S_0 = 0;
    %% CSA (4a)
    str_csa = 'sqN';
    MODE_CSA = 1;
    E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    C_handle = @(N,MU) 1/sqrt(N);
    D_handle = @(N,MU) 1/C_handle(N,MU);
    %% CSA (4b)
    % str_csa = 'linN';
    % MODE_CSA = 1;
    % E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    % C_handle = @(N,MU) 1/N;
    % D_handle = @(N,MU) 1/C_handle(N,MU);
    %% CSA (4c)
    % str_csa = 'han';
    % MODE_CSA = 2;
    % E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    % C_handle = @(N,MU) (MU+2)/(MU+N+5);
    % D_handle = @(N,MU) 1 + C_handle(N,MU) + 2 * max(0, sqrt((MU-1)/(N+1))-1);
elseif strcmpi(MODE, 'SA')
    tau_handle = @(N) 1/sqrt(2*N); %1/sqrt(N)*sqrt((1-0.5^2)/(1-1/(2*0.8*sqrt(2*N)))) 
end


%% Y via R
R_0 = 1e9;
SIGMA_STOP = 1e-15;
R_STOP = 1e-12;
F_STOP = nan; %1e-3;
% SIGMA_0: initilaized as fct of MU

%% Initialize

flag_success = nan*zeros(TRIALS, 1);
feval_t = zeros(TRIALS, 1);
r_g_t = zeros(G, TRIALS); 
f_g_t = zeros(G, TRIALS); 
sigma_g_t = zeros(G, TRIALS); 
signorm_g_t = zeros(G, TRIALS); 
%y_g_t = zeros(N, G, TRIALS); 


for j=1:length(LIST_SIGMA_RESCALE)
    SIGMA_RESCALE = LIST_SIGMA_RESCALE(j);

    %% Main Loop 
    FIT = Fitness('Sphere',N,[]);
    Y_0 = R_0/sqrt(N)*ones(N,1);
    SIGNORM_ZERO = sqrt( sqrt(N*(8*e_vartheta_a_b(THETA,1,0)^2*MU_SWEEP(1)^2 + N)) - N);
    SIGMA_0 = SIGNORM_ZERO * norm(Y_0) / N;  
    
    rng(1);
    for i=1:TRIALS
        fprintf('  t = %i \n', i)
        if strcmpi(MODE, 'SA')
            % TAU = tau_handle(N);
            % [y, f_g, r_g, sigma_g, g_end, ~] = muComLam_sSA(FIT, N, MU, LAM, Y_0, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, G, verbose);
            error('sSA undefined mu sweep.');
        elseif strcmpi(MODE, 'CSA')
            E_chi = E_chi_handle(N);
            [y, f_g, r_g, sigma_g, g_end, ~] = csa_es(FIT, N, MU_SWEEP, THETA, Y_0, SIGMA_0, S_0, C_handle, D_handle, E_chi, ...
                SIGMA_STOP, R_STOP, F_STOP, G, verbose, MODE_CSA, SIGMA_RESCALE,j);
        end
    
        f_g_t(:,i) = f_g;
        sigma_g_t(:,i) = sigma_g;
        r_g_t(:,i) = r_g;
        signorm = sigma_g./r_g*N;
        signorm_g_t(:,i) = signorm;
    
    end
end % SIGMA_RESCALE

xticks([0:500:2000]);
yticks(10.^(-12:3:30));
myfigsize(gcf,6,3.5,8,8);
str = ['N',num2str(N),'_',str_csa]; 
% savefig(gcf,[str, '.fig']); 
% exportgraphics(gcf,[str, '.pdf'], 'ContentType','vector');



