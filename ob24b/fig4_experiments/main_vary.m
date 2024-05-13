clear
addpath('../code');

%% SET log-normal (0) or normal (1) mutations
SET_metaEP = 0;


%% General Settings
verbose = 0;
MODE = 'SA';
tau_fct = @(N) 1/sqrt(2*N);
TRIALS = 10;

%% DOE1
MU = 500;
LAM = 2000;
DOE = 1;
N_LIST = [20];
A_LIST = [1,2,5,10,20,50,100];
ALPHA_LIST = [2*pi, 4*pi, 6*pi, 8*pi, 10*pi, 12*pi];
FEVAL_STOP = 12e7;
%EP needs on avg. 12e6 feval for local conv. on N=20,A=100,Alpha=12pi => fac. 10

%% DOE2
% MU = 1000;
% LAM = 2000;
% DOE = 2;
% N_LIST = [100];
% A_LIST = [1,2,5,10,20,50]; %[1,2,5,10,20,50,100];
% ALPHA_LIST = [pi, 2*pi, 4*pi, 6*pi, 8*pi]; %[pi, 2*pi, 4*pi, 6*pi, 8*pi];
% FEVAL_STOP = 20e7; 
%EP needs on avg. 20e6 feval for local conv. on N=100,A=50,Alpha=8pi => fac. 10
%EP needs on avg. 60e6 feval for local conv. on N=100,A=100,Alpha=10pi

N_dim = length(N_LIST);
N_A = length(A_LIST);
N_alpha = length(ALPHA_LIST);

%% RNG
rng(1); 

%% Stop
SIGMA_STOP = 1e-4;
SIGMA_SS_TOL = 0.1;
R_STOP = 1e-3;
F_STOP = nan;

%% Main Loop 
tic
G = FEVAL_STOP/LAM+1;
N_VAR = N_dim*N_A*N_alpha;
init = zeros(N_dim,N_A,N_alpha);
[dat_N,dat_A,dat_ALPHA,dat_feval,dat_suc,dat_f_delta,dat_sign_ss] = ...
    deal(init,init,init,init,init,init,init);

c = 0;
for n=1:N_dim
    %% N-loop
    N = N_LIST(n);
    TAU = tau_fct(N);
    %fprintf('\t n = %i: \n', N);
    Y_HAT = zeros(N,1);
    
    for t=1:N_A
    A = A_LIST(t);
        
    for v=1:N_alpha

        %fprintf('\t theta = %d: \n', VARTHETA);
        ALPHA = ALPHA_LIST(v);
        Y_0 = ceil(ALPHA*A/2)*ones(N,1); 
        FIT = @(x) sum(x.^2 + A*(1 - cos(ALPHA * x)), 1);
          
        %% VARIATION
        c = c+1;
        fprintf('\t var = %i/%i: \n', c, N_VAR);
        
        %% MU-loop
        %LAM = LAM_LIST(m);
        %fprintf('\t mu = %i: \n', MU);
        %MU = round(LAM*VARTHETA);
        
        %% Inner loop 
        e_10 = e_vartheta_a_b(MU/LAM,1,0);
        x_2ndZero = sqrt( sqrt(N*(8*e_10^2*MU^2 + N)) - N);
        SIGMA_0 = x_2ndZero * norm(Y_0) / N;   
        sigma_eps = sqrt(N/2)*A;
        R_infty = sqrt(sigma_eps*N/(4*e_10*MU));
        sign_ss = sqrt(MU);
        sigma0 = sqrt(sigma_eps/(4*e_10*N)); 
    
        for i = 1:TRIALS
            fprintf('\t\t t = %i:', i);

            %% Single Run
            if strcmpi(MODE, 'SA')
                [y, f_g, r_g, sigma_g, gen, feval] = muComLam_metaEP_noFclass(SET_metaEP, FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, FEVAL_STOP, G, verbose);
            end

            %% Data
            dat_N(n,t,v,i) = N;
            dat_A(n,t,v,i) = A;
            dat_ALPHA(n,t,v,i) = ALPHA;
            dat_feval(n,t,v,i) = feval;

            if r_g(gen) < R_STOP
                dat_suc(n,t,v,i) = 1;
            elseif sigma_g(gen) < SIGMA_STOP
                dat_suc(n,t,v,i) = 0;
            elseif feval == FEVAL_STOP
                if abs(mean(sigma_g(gen-1000:gen)-sigma0))/sigma0 < SIGMA_SS_TOL %abs((mean(sigma_g(4000:5000))-mean(sigma_g(3000:4000)))/mean(sigma_g(3000:4000)))<0.05
                    dat_suc(n,t,v,i) = -1;
                else
                    warning('FEVAL_STOP without SteadyState');
                    dat_suc(n,t,v,i) = nan;
                end
            end
            g0 = round(gen/5);
            dat_sign_ss(n,t,v,i) = mean(sigma_g(g0:gen)./r_g(g0:gen)*N);

        end % trials 
    
    end %vartheta
    end %N
    
    %% save data
    %save('data_wsp.mat', 'data');
    
end %tau
toc
if SET_metaEP==1, SA_MODE = 'EP'; else, SA_MODE = 'LN'; end


%% CHECK ERT

%% SAME COLORBAR SCALE FOR ALL PLOTS (CONTOUR)

%% END: main loop

ERT = nan*zeros(N_dim,N_A,N_alpha);
ECDF = cell(N_dim,N_A,N_alpha);
P_S = nan*zeros(N_dim,N_A,N_alpha);
for n=1:N_dim
    for t=1:N_A
        for v=1:N_alpha
            %for m=1:N_lam
                ids_suc = find([dat_suc(n,t,v,:)]==1);
                ids_unsuc = find([dat_suc(n,t,v,:)]==0);
                num_suc = length(ids_suc);
                num_unsuc = length(ids_unsuc);

                %% ECDF
                % sort all feval values to generate increasing ECDF
                [val_sorted_feval,ids_sorted_feval] = sort([dat_feval(n,t,v,ids_suc)]);
                % calculate ECDF (increasing for each value)
                if ~isempty(ids_sorted_feval)
                    ECDF{n,t,v} = (1:length(ids_sorted_feval))/TRIALS;
                end   
                P_S(n,t,v) = num_suc/TRIALS;

                %% ERT
                if ~isempty(ids_suc)
                    ERT(n,t,v) = sum([dat_feval(n,t,v,:)])/num_suc;
                else
                    ERT(n,t,v) = nan;
                end
                
                %plot(val_sorted_feval, ECDF{t,n,v}, '-', 'DisplayName', [var_str, '=', num2str(VAR_LIST(j))]);
  
            %end
        end
    end
end

str_save = ['DOE',num2str(DOE),'_SA', SA_MODE, '_TR', num2str(TRIALS)];
save([str_save, '.mat']);
