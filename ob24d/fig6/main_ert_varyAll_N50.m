clear
addpath('../code');

verbose = 0;
G = 1e5+1;
MODE = 'SA';

%% Variation
TRIALS = 200;
N_LIST = [50];

%% LARGE DOE
LAM_LIST = [200:200:2000]; %[100:100:1000];
TAU_LIST = [0.05:0.02:0.17]; %0.06:0.02:0.20
% TAU_LIST = [eval('1/sqrt(8*N_LIST)'),eval('1/sqrt(6*N_LIST)'),eval('1/sqrt(4*N_LIST)'),eval('1/sqrt(3*N_LIST)'),eval('1/sqrt(2*N_LIST)'),eval('1/sqrt(1.5*N_LIST)'),eval('1/sqrt(1.1*N_LIST)'),eval('1/sqrt(N_LIST)')]; %{'1/sqrt(2*N)','1/sqrt(N)'};
THETA_LIST = [0.1:0.1:0.7]; %[0.1:0.1:0.7];

%% MEDIUM DOE
% LAM_LIST = [200:100:2000];
% TAU_LIST = [eval('1/sqrt(8*N_LIST)'),eval('1/sqrt(4*N_LIST)'),eval('1/sqrt(2*N_LIST)')];
% THETA_LIST = [0.2:0.1:0.6];
%% SINGLE VAR
% LAM_LIST = [1000];
% TAU_LIST = [eval('1/sqrt(8*N_LIST)')];
% THETA_LIST = [0.3];

N_dim = length(N_LIST);
N_lam = length(LAM_LIST);
N_tau = length(TAU_LIST);
N_theta = length(THETA_LIST);

%% RNG
rng(1); 

%% Stop
SIGMA_STOP = 1e-3;
R_STOP = nan;
F_STOP = 1e-1;
FEVAL_STOP = 1e7;

 %% Rastrigin Fitness
A = 10;
ALPHA = 2*pi;
FIT = @(x) sum(x.^2 + A*(1 - cos(ALPHA * x)), 1);

%% Main Loop 
tic
N_VAR = N_dim*N_tau*N_theta*N_lam;
init = zeros(N_dim,N_tau,N_theta,N_lam);
[dat_N,dat_MU,dat_TAU,dat_THETA,dat_feval,dat_suc,dat_f_delta,dat_sign_ss] = ...
    deal(init,init,init,init,init,init,init,init);

c = 0;
for n=1:N_dim
    %% N-loop
    N = N_LIST(n);
    %fprintf('\t n = %i: \n', N);
    Y_0 = 30*ones(N,1); 
    Y_HAT = zeros(N,1);
    
    for t=1:N_tau
    TAU = TAU_LIST(t);
    %fprintf(['\t tau = ',TAU, '\n']);
        
     for v=1:N_theta
        VARTHETA = THETA_LIST(v);
        %fprintf('\t theta = %d: \n', VARTHETA);
    
    for m=1:N_lam
        
        %% VARIATION
        c = c+1;
        fprintf('\t var = %i/%i: \n', c, N_VAR);
        
        %% MU-loop
        LAM = LAM_LIST(m);
        %fprintf('\t mu = %i: \n', MU);
        MU = round(LAM*VARTHETA);
        e_10 = e_mu_lam_a_b_v2(MU,LAM,1,0);
        
        %% Inner loop 
        x_2ndZero = sqrt( sqrt(N*(8*e_10^2*MU^2 + N)) - N);
        SIGMA_0 = x_2ndZero * norm(Y_0) / N;        
    
        for i = 1:TRIALS
            fprintf('\t\t t = %i:', i);

            %% Single Run
            if strcmpi(MODE, 'SA')
                [y, f_g, r_g, sigma_g, gen, feval] = muComLam_sSA(FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, FEVAL_STOP, G, verbose);
            end

            %% Data
            dat_N(n,t,v,m,i) = N;
            dat_MU(n,t,v,m,i) = MU;  
            dat_TAU(n,t,v,m,i) = TAU;
            dat_THETA(n,t,v,m,i) = VARTHETA;
            dat_feval(n,t,v,m,i) = feval;

            if f_g(gen) < F_STOP
                dat_suc(n,t,v,m,i) = 1;
                dat_f_delta(n,t,v,m,i) = 0;
            else
                dat_suc(n,t,v,m,i) = 0;
                dat_f_delta(n,t,v,m,i) = log10(f_g(gen)/F_STOP);
            end
            g0 = round(gen/5);
            dat_sign_ss(n,t,v,m,i) = mean(sigma_g(g0:gen)./r_g(g0:gen)*N);

        end % trials 
    
    end %mu
    end %vartheta
    end %N
    
    %% save data
    %save('data_wsp.mat', 'data');
    
end %tau
toc
save('wsp.mat');

%% CHECK ERT

%% SAME COLORBAR SCALE FOR ALL PLOTS (CONTOUR)

%% END: main loop

ERT = nan*zeros(N_dim,N_tau,N_theta,N_lam);
ECDF = cell(N_dim,N_tau,N_theta,N_lam);
P_S = nan*zeros(N_dim,N_tau,N_theta,N_lam);
for n=1:N_dim
    for t=1:N_tau
        for v=1:N_theta
            for m=1:N_lam
                ids_suc = find([dat_suc(n,t,v,m,:)]==1);
                ids_unsuc = find([dat_suc(n,t,v,m,:)]==0);
                num_suc = length(ids_suc);
                num_unsuc = length(ids_unsuc);

                %% ECDF
                % sort all feval values to generate increasing ECDF
                [val_sorted_feval,ids_sorted_feval] = sort([dat_feval(n,t,v,m,ids_suc)]);
                % calculate ECDF (increasing for each value)
                if ~isempty(ids_sorted_feval)
                    ECDF{n,t,v,m} = (1:length(ids_sorted_feval))/TRIALS;
                end   
                P_S(n,t,v,m) = num_suc/TRIALS;

                %% ERT
                if ~isempty(ids_suc)
                    ERT(n,t,v,m) = sum([dat_feval(n,t,v,m,:)])/num_suc;
                else
                    ERT(n,t,v,m) = nan;
                end
                
                %plot(val_sorted_feval, ECDF{t,n,v,m}, '-', 'DisplayName', [var_str, '=', num2str(VAR_LIST(j))]);
  
            end
        end
    end
end