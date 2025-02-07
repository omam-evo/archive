clear


THETA = 0.5;
LIST_CONFIG = [0];

%% Approx 1a: highest quality, N-dep shows deviations at small N
TAYLOR_MODE = 0; %0,1,2,3
SET_RATIO_APPROX_sA = 0; % 0,1
PHI_MODE = 'full';  %'full','medium','large', 'arn','large_v2'
MODE_RES_SIGNZERO = 'full'; %'full','large'

%% Approx 1b: ratio 1 for sA
% TAYLOR_MODE = 0; %0,1,2,3
% SET_RATIO_APPROX_sA = 1; % 0,1
% PHI_MODE = 'full';  %'full','medium','large', 'arn','large_v2'
% MODE_RES_SIGNZERO = 'full'; %'full','large'

%% Approx 2a: ratio 1 for sA
% TAYLOR_MODE = 2; %0,1,2,3
% SET_RATIO_APPROX_sA = 1; % 0,1
% PHI_MODE = 'large';  %'full','medium','large', 'arn','large_v2'
% MODE_RES_SIGNZERO = 'large'; %'full','large'

%% Approx 2b: ratio 1 for sA
% TAYLOR_MODE = 3; %0,1,2,3
% SET_RATIO_APPROX_sA = 1; % 0,1
% PHI_MODE = 'large';  %'full','medium','large', 'arn','large_v2'
% MODE_RES_SIGNZERO = 'large'; %'full','large'

%% Approx 3: ratio 1 for sA
% TAYLOR_MODE = 0; %0,1,2,3
% SET_RATIO_APPROX_sA = 1; % 0,1
% PHI_MODE = 'large';  %'full','medium','large', 'arn','large_v2'
% MODE_RES_SIGNZERO = 'large'; %'full','large'

%% Approx to solve Eq.
% TAYLOR_MODE = 2; %0,1,2
% SET_RATIO = 1; % 0,1
% PHI_MODE = 'large';  %'full','medium','large', 'arn'

%% CHOOSE CSA
C_handle = @(MU,N) 1/sqrt(N);
D_handle = @(MU,N) 1/C_handle(MU,N);

% C_handle = @(MU,N) 1/N;
% D_handle = @(MU,N) 1/C_handle(MU,N);
 
% C_handle = @(MU,N) (MU+2)/(MU+N+5);
% D_handle = @(MU,N) 1 + 1/C_handle(MU,N) + 2*max(0, sqrt((MU-1)/(N+1))-1)/C_handle(MU,N);

for k=1:length(LIST_CONFIG)
    CONFIG = LIST_CONFIG(k);
    fprintf('CONFIG=%i\n',CONFIG)

    if CONFIG == 0 %% N VARIATION
        N_LIST = [100];
        MU_LIST = 1000*ones(1,length(N_LIST));
        TRIAL_LIST = [1]; %1*ones(1,length(N_LIST));
        RES_SIGN_OPT = [nan];
        RES_SIGN_ZERO = [nan];
        MODE_MU_N = 1;
        configstr = 'test';

    elseif CONFIG == 1 %% MU var, N small
        MU_LIST = [10,30,100,300,1000,3000,10000];
        N_LIST = 10*ones(1,length(MU_LIST));
        RES_SIGN_OPT = nan*ones(1,length(N_LIST)); %[4.0006    5.4007    8.2508   11.4008   15.6008   18.0009];
        RES_SIGN_ZERO = [8.9401   16.3441   30.4921   53.0671   97.2401  168.6601 307.3343];
        TRIAL_LIST = [100,100,100,100,40,40,40];
        MODE_MU_N = 1;
        configstr = 'conf1';

    elseif CONFIG == 2
        MU_LIST = [10,30,100,300,1000,3000,10000];
        N_LIST = 100*ones(1,length(MU_LIST));
        RES_SIGN_OPT = nan*ones(1,length(MU_LIST)); %[6.30055	11.2006	18.55065 27.9007 42.330751 44.25085];
        RES_SIGN_ZERO = nan*ones(1,length(MU_LIST)); %[12.16613	24.92011	47.9651	84.35109	154.70009	268.45009];     
        TRIAL_LIST = [100,100,100,100,40,40,40];
        MODE_MU_N = 1;
        configstr = 'conf2';

    elseif CONFIG == 3 %% N VARIATION
        N_LIST = [10,20,50,100,200,500,1000]; %,500,1000
        MU_LIST = 1000*ones(1,length(N_LIST));
        RES_SIGN_OPT = nan*ones(1,length(N_LIST)); %[15.6008500000000	17.8508500000000	29.2008000000000	34.0008000000000	50.2507500000000	62.5007500000000	88.5007000000000];
        RES_SIGN_ZERO = [96.8241	109.5991	132.5681	nan nan nan nan];
        TRIAL_LIST = [100,100,100,50,25,10,5]; %1*ones(1,length(N_LIST));
        MODE_MU_N = 2;
        configstr = 'conf3';

    elseif CONFIG == 4
        N_LIST = [10,20,50,100,200,500,1000];
        MU_LIST = 2*N_LIST;
        RES_SIGN_OPT = nan*ones(1,length(MU_LIST));
        RES_SIGN_ZERO = [13.122	21.393	41.295 nan nan nan nan];
        TRIAL_LIST = [100,100,100,50,25,10,5];
        MODE_MU_N = 2;
        configstr = 'conf4';
    end

for i=1:length(MU_LIST)
    MU = MU_LIST(i);
    N = N_LIST(i);
    LAM = MU/THETA;

    G = 1000;
    C = C_handle(MU,N);
    D = D_handle(MU,N);
    E_chi = sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
    sfac = sqrt(C*(2-C)*MU);
    sph = Sphere;
    
    phi_g = nan*zeros(G,1);
    sa = phi_g;
    ss = phi_g;
    sign = phi_g;
    
    sa(1) = 0;
    ss(1) = 0;
    
  if strcmpi(PHI_MODE,'full') || strcmpi(PHI_MODE,'arn')
    e_10 = e_mu_lam_a_b_v2(MU,LAM,1,0);
    e_11 = e_mu_lam_a_b_v2(MU,LAM,1,1);
    e_20 = e_mu_lam_a_b_v2(MU,LAM,2,0);
  else %'medium', 'large'
    e_10 = e_vartheta_a_b(MU/LAM,1,0);
    e_11 = e_vartheta_a_b(MU/LAM,1,1);
    e_20 = e_vartheta_a_b(MU/LAM,2,0);
  end

    sign(1) = (8*N)^(1/4)*sqrt(e_10*MU);
    if strcmpi(PHI_MODE,'full') || strcmpi(PHI_MODE,'arn') || strcmpi(PHI_MODE,'medium') 
        E_za = @(x) e_10/sqrt(1+x^2/(2*N));
    elseif strcmpi(PHI_MODE,'large') || strcmpi(PHI_MODE,'large_v2')
        E_za = @(x) e_10*sqrt(2*N)/x;
    end
    if strcmpi(PHI_MODE,'full') || strcmpi(PHI_MODE,'arn')
        E_zz = @(x) N/MU*(1 + 1/N*(e_11+(MU-1)*e_20)./(1+x.^2/(2*N)) - (N-1)/N^2.*x*e_10./sqrt(1+x.^2/(2*N)));
    elseif strcmpi(PHI_MODE,'medium') || strcmpi(PHI_MODE,'large_v2') 
        E_zz = @(x) N/MU*(1 + 2*(e_11+MU*e_20)/x.^2 - e_10*sqrt(2/N));
    elseif strcmpi(PHI_MODE,'large')
        E_zz = @(x) N/MU;
    end

    %% NUMERIC SOLVING
    % ratio = @(x) 1 - phi(x, MU, N, {PHI_MODE}, e_10, e_11, e_20)/N;
    % rhs = @(x) MU*zz + 2*(1-C)*MU*(E_za(x).^2-x/N*E_za(x)*zz)/(C+(1-C)*x/N*E_za(x));
    % 
    % f0 = @(x) E_chi^2*(1+D*log(ratio(x))).^2 - rhs(x);
    % f1 = @(x) E_chi^2*(1-phi(x, MU, N, {PHI_MODE}, e_10, e_11, e_20)/N).^2 - rhs(x);
    % f2 = @(x) N*(1-2*phi(x, MU, N, {PHI_MODE}, e_10, e_11, e_20)/N) - rhs(x);
    
    
    for g=1:G
        za = E_za(sign(g));
        zz = E_zz(sign(g));
        ratio_Rgp_Rg = 1 - phi(sign(g), MU, N, {PHI_MODE}, e_10, e_11, e_20)/N;
       
        %% SET RATIO
        if SET_RATIO_APPROX_sA == 0
            ratio_sA = ratio_Rgp_Rg;  
        else
            ratio_sA = 1;
        end
        %% Eq. (6.30): iteration of sA in expectation, without covariance
        sa(g+1) = 1/ratio_sA * ((1-C)*sa(g) - (1-C)*sign(g)/N*sa(g)*za + sfac*za - sfac*sign(g)/N*zz );

    
        %% Eq. (6.25)
        ss(g+1) = (1-C)^2*ss(g) + 2*(1-C)*sfac*sa(g)*za + sfac^2*zz;
    
        %% full
        if TAYLOR_MODE==0
            sign(g+1) = sign(g)*1/ratio_Rgp_Rg*exp(1/D*(sqrt(ss(g+1))/E_chi - 1));
        elseif TAYLOR_MODE==1
            sign(g+1) = sign(g) * ...
                (1+phi(sign(g), MU, N, {PHI_MODE}, e_10, e_11, e_20)/N) * ...
                (1+(sqrt(ss(g+1))-E_chi)/(D*E_chi));
        elseif TAYLOR_MODE==2
            sign(g+1) = sign(g) * ... 
                (1 + phi(sign(g), MU, N, {PHI_MODE}, e_10, e_11, e_20)/N + ... 
                (sqrt(ss(g+1))-E_chi)/(D*E_chi));
        elseif TAYLOR_MODE==3
            sign(g+1) = sign(g) * ... 
                (1 + phi(sign(g), MU, N, {PHI_MODE}, e_10, e_11, e_20)/N + (ss(g+1)-N)/(2*D*N));
        end

        assert(isreal(sqrt(ss(g+1))));
    end
    
    if length(MU_LIST)==1
        figure; hold on; legend;
            plot(sa, 'r-','DisplayName','$s_A$')
            plot(ss, 'b-','DisplayName','$||\mathbf{s}||^2$')
            plot(sign, 'k-','DisplayName','$\sigma^*$');
            xlabel('$g$');
            ylabel('Dynamics');
            xscale('log');
    end
    assert(sign(end)-sign(end-1)<1e-3);
    
    %% normalize using sigma*_0
    if isnan(RES_SIGN_ZERO(i)) && strcmpi(PHI_MODE,'full')
        [RES_SIGN_OPT(i), RES_SIGN_ZERO(i)] = ...
                sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU,LAM,1,0), MU, LAM, N, {'full'});      
    elseif strcmpi(PHI_MODE,'large') %% normalize using sigma*_phi0 (analytic)
        RES_SIGN_OPT(i) = 0;
        RES_SIGN_ZERO(i) = (8*N)^(1/4)*sqrt(e_10*MU);
    end
 
    sign_iter(i) = sign(end);
    signzero_approx = (8*N)^(1/4)*sqrt(e_10*MU);
    
    %% NUMERIC SOLVING
    % sign_sol_f0 = fzero(f0,signzero_approx);
    % sign_sol_f1 = fzero(f1,signzero_approx);
    % sign_sol_f2 = fzero(f2,signzero_approx);
    sign_sol_eq(i) = csa_predict_sign(MU,LAM,N,C,D,0);  %uses approx (6.41,6.42,6.43)
    gamma(i) = csa_get_gamma(THETA,C,D,N);              %uses approx of sign_sol_eq, and uses approx. gamma*signzero
    % gamma_f0(i) = sign_sol_f0/RES_SIGN_ZERO(i);
    % gamma_f1(i) = sign_sol_f1/RES_SIGN_ZERO(i);
    % gamma_f2(i) = sign_sol_f2/RES_SIGN_ZERO(i);
    % gamma_eq(i) = sign_sol_eq/RES_SIGN_ZERO(i);
end
% if MODE_MU_N == 1
%     x = MU_LIST;
%     %xstr = '$\mu$';
% else
%     x = N_LIST;
%     %xstr = '$N$';
% end
% if k==1
%     figure;
% else
%     hold on;
% end
% 
%     plot(x,sign_iter./RES_SIGN_ZERO, 'DisplayName',['iter: config:',num2str(CONFIG)]);
%     %plot(x,sign_sol_eq./RES_SIGN_ZERO, 'k--', 'DisplayName','solve sign');
%     %plot(x,gamma, 'k-.','DisplayName','solve gamma');
%     xlabel('Either MU or N');
%     ylabel('GAMMA')
%     set(gca,'Xscale','log')

end

% x = linspace(0,1.1*signzero_approx,1001);
% phi_arn = sph.phi_arn(x, MU, N, e_10, e_11, e_20);
% phi = sph.phi(x, e_10, MU, LAM, N, {'full'});
% figure; hold on;
%     plot(x, phi)
%     plot(x, phi_arn)

% sign_full = signzero;
% [~,sign_medium] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU,LAM,1,0), MU, LAM, N, {'medium'});
% sign_approx = (8*N)^(1/4)*sqrt(e_10*MU);