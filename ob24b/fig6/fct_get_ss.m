function [sign_ss_res, ss_path, C_handle, D_handle, r_g, sigma_g] = fct_get_ss(THETA, MODE, TYPE,N_LIST,MU_LIST,TRIAL_LIST,G_MAX,SET_META_EP)

THETA = 0.5;

verbose = 0;

%% Misc
if strcmpi(MODE, 'CSA')
    S_0 = 0;
    if strcmpi(TYPE, 'AR')
        MODE_CSA = 0;
        E_chi_handle = @(N) N;
        C_handle = @ (N,MU) 1/sqrt(N);
        D_handle = @(N,MU) 1/C_handle(N,MU);
    elseif strcmpi(TYPE, 'HanV1') %% CONFIG A1: Hansen derandomized 1/sqrt
        MODE_CSA = 1;
        E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
        C_handle = @(N,MU) 1/sqrt(N);
        D_handle = @(N,MU) 1/C_handle(N,MU);
    elseif strcmpi(TYPE, 'HanV1_linN') %% CONFIG A2: Hansen derandomized
        MODE_CSA = 1;
        E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
        C_handle = @(N,MU) 1/N;
        D_handle = @(N,MU) 1/C_handle(N,MU);
    elseif strcmpi(TYPE, 'HanV2') %% CONFIG B: Hansen tutorial default
        MODE_CSA = 2;
        E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
        C_handle = @(N,MU) (MU+2)/(MU+N+5);
        D_handle = @(N,MU) 1 + C_handle(N,MU) + 2 * max(0, sqrt((MU-1)/(N+1))-1);
    elseif strcmpi(TYPE, 'HanV2_mod') %% Hansen tutorial modified (mimic derandomized 1/sqrt)
        MODE_CSA = 2;
        E_chi_handle = @(N) sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
        C_handle = @(N,MU) 1/sqrt(N);
        D_handle = @(N,MU) 1 + C_handle(N,MU);   
    else
        error('CSA TYPE not defined.')
    end
elseif strcmpi(MODE, 'SA')
    C_handle = @(x,y) nan;
    D_handle = @(x,y) nan;
    s_g = nan;
    if strcmpi(TYPE, 'N')
        tau_handle = @(N) 1/sqrt(N);
    elseif strcmpi(TYPE, '2N')
        tau_handle = @(N) 1/sqrt(2*N);
    elseif strcmpi(TYPE, '8N')
        tau_handle = @(N) 1/sqrt(8*N);
    end
end

%% Y via R
SIGMA_STOP = 1e-10;
R_STOP = nan;
F_STOP = 1e-6; %1e-3;
% SIGMA_0: initilaized as fct of MU

%% Initialize
% sph = Sphere;
% mode = 'full';
sign_ss_res = nan*zeros(1,length(MU_LIST));
% signzero_approx = nan*sign_ss_res;
% sign_phi0_full = nan*sign_ss_res;
% sign_opt = nan*sign_ss_res;
gen = zeros(length(MU_LIST),1);

ss_path = sign_ss_res;  % ||s||^2
sign_pred = sign_ss_res;    % analytic formula for sigma^*
sign_ratio = sign_ss_res;   % ratio of sigma^*_ss to one-gen second zero sigma^*_0
signzero_analytic = sign_ss_res; % signzero from varphi^*

%% Main Loop 
tic
for m=1:length(MU_LIST)
    
    MU = MU_LIST(m);
    N = N_LIST(m);
    LAM = round(MU/THETA);
    fprintf('mu = %i, N = %i \n', MU, N)

    G = G_MAX(N);

    %% TRIALS
    TRIALS = TRIAL_LIST(m);
    flag_success = nan*zeros(TRIALS, 1);
    feval_t = zeros(TRIALS, 1);
    r_g_t = zeros(G, TRIALS); 
    f_g_t = zeros(G, TRIALS); 
    sigma_g_t = zeros(G, TRIALS); 
    signorm_g_t = zeros(G, TRIALS); 
    sign_ss_per_trial = nan*zeros(length(MU_LIST),TRIALS);

    
    %% General
    FIT = Fitness('Sphere',N,[]);
    Y_0 = 1*ones(N,1);
    SIGNORM_ZERO = (8*N)^(1/4)*sqrt(e_vartheta_a_b(MU/LAM,1,0)*MU);
    SIGMA_0 = SIGNORM_ZERO * norm(Y_0) / N;  
    signzero_analytic(m) = SIGNORM_ZERO;

    rng(1);
    FLAG_SIGMA_STOP = 0;
    for i=1:TRIALS
        fprintf('  t = %i \n', i)
        if strcmpi(MODE, 'SA')
            TAU = tau_handle(N);
            [y, f_g, r_g, sigma_g, g_end, ~] = muComLam_sSA(SET_META_EP,FIT, N, MU, LAM, Y_0, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, G, verbose);
        elseif strcmpi(MODE, 'CSA')
            E_chi = E_chi_handle(N);
            C = C_handle(N,MU);
            D = D_handle(N,MU);
            [y, f_g, r_g, sigma_g, g_end, s_g] = csa_es(FIT, N, MU, LAM, Y_0, SIGMA_0, S_0, C, D, E_chi, SIGMA_STOP, R_STOP, F_STOP, G, 0, MODE_CSA);
        end

        f_g_t(:,i) = f_g;
        sigma_g_t(:,i) = sigma_g;
        r_g_t(:,i) = r_g;
        signorm = sigma_g./r_g*N;
        signorm_g_t(:,i) = signorm;
        gen(m) = gen(m) + g_end;
        % sign_ss_per_trial(m,i) = mean(signorm(100:end), 'omitnan'); 

        % if sigma_g(g_end)<SIGMA_STOP
        %     signorm_g_t(:,i) = nan*signorm;
        % end
        if sigma_g(g_end)<SIGMA_STOP
            FLAG_SIGMA_STOP = FLAG_SIGMA_STOP+1;
            % break;
        end

    end

    %% Experiment
    if FLAG_SIGMA_STOP==0
        f_g_med = median(f_g_t, 2,'includenan');
        sign_g_med = median(signorm_g_t, 2,'includenan');
        signorm_g = sign_g_med(~isnan(sign_g_med)); % select all not-nan values 
        G_0 = round(length(signorm_g)/2);
        sign_ss_res(m) = median(signorm_g(G_0:end));
    
        ss = vecnorm(s_g,2,1).^2;
        ss_path(m) = median(ss(G_0:end), 'omitnan');
    else
        sign_ss_res(m) = nan;
        ss_path(m) = nan;
    end
    %% DEBUG 1
    % figure; hold on;
    %     ids=[1,2,3,56,58,65];
    %     for ii=ids
    %         plot(r_g_t(:,ii), 'k-')
    %         plot(sigma_g_t(:,ii), 'r--')
    %     end


    %% DEBUG 1
    % figure; plot(signorm_g)

    %% DEBUG SIGN-DYN
    % figure; hold on; 
    %     plot(signorm_g_t); 
    %     plot(signorm_g, 'k-', 'LineWidth',2)
    % a = 1;

    %% What is this?
    % if length(LIST_MU)==1
    %     figure; hold on;
    %         plot(ss, 'k-');
    %         plot([0,g_end],[sign_ss_res(m),sign_ss_res(m)],'r-','DisplayName','meas');
    %         plot([0,g_end],[res,res],'b-', 'DisplayName', 'pred');
    % 
    % end

    % sign_pred(m) = csa_predict_sign(MU,LAM,N,C,D);

end % MU
toc

end % function
