function main_iter(MU, LAM,A, N, ALPHA, SIGMA_NORM, R_STOP, F_STOP, N_G, SEED, ITER_METHOD )
%% PARAMS
% MU = 100;
% LAM = 200;
% A = 1;
% N = 100;
% ALPHA = 2*pi;
% sph = Sphere;
% SIGMA_NORM = 5; %1*sph.signorm_opt_sph(C_MU_LAM, MU, LAM, N, {'medium'});
% F_STOP = nan; % 2e-6;
% R_STOP = 1e-6;  %1e-4;
% N_G = 1000;
% SEED = 1;


%% Iteration method
% ITER_METHOD = 'Y'; % 'SIM' using simulations, 'R' via phi^(II)(R),'Y' via phi^(II),'PHI1' via phi^(1)
C_MU_LAM = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
c_vartheta = e_vartheta_a_b(MU/LAM, 1, 0);
e_11 = 0*e_mu_lam_a_b_v2(MU, LAM, 1, 1); % Adjust to include terms for phi^(II)(y_i)
e_20 = 0*e_mu_lam_a_b_v2(MU, LAM, 2, 0); % Adjust to include terms for phi^(II)(y_i)

%% INIT
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1);
FIT_R = @(R) R^2 + N*A*(1 - exp(-0.5*(ALPHA*R/sqrt(N))^2));
ras_tools = RasTools;
Phi = Phi_v2;
[order_y_min, order_y_max] = ras_tools.get_all_extrema(ALPHA, A, N);

%% Position
a = 20;   % R > 7*sqrt(N) for A=10 necessary
R_0 = a*sqrt(N);
rng(SEED); v = randn(N,1); 
y = v/norm(v)*R_0;

%% Data
y_g = nan*zeros(N, N_G+1);
k_g = nan*zeros(N, N_G+1);
d_g = nan*zeros(N, N_G+1);
R_g = nan*zeros(N_G+1, 1);
F_g = nan*zeros(N_G+1, 1);
sigma_g = nan*zeros(N_G+1, 1);
D_Q_g = nan*zeros(N, N_G+1);

if strcmp(ITER_METHOD, 'R') || strcmp(ITER_METHOD, 'TEST')
    R = R_0;
    phi_g = nan*zeros(1, N_G+1);
else
    phi_g = nan*zeros(N, N_G+1);
end

for i=1:N_G+1
    
    if strcmp(ITER_METHOD, 'R') || strcmp(ITER_METHOD, 'TEST')
        %% R iteration
        R_g(i) = R;
        F_g(i) = FIT_R(R);
        sigma = SIGMA_NORM*R/N;
        sigma_g(i) = sigma; 
    else %'SIM', 'Y', 'PHI1'
        %% y_i iteration
        R = norm(y);
        sigma = SIGMA_NORM*R/N;

        R_g(i) = R;
        F_g(i) = FIT(y);
        sigma_g(i) = sigma;
        y_g(:,i) = y;
    end

    fprintf('\t Iter %i, R(g)=%d, F(g)=%d, sigma(g)=%d \n', i, R_g(i), F_g(i), sigma_g(i));

    if strcmp(ITER_METHOD, 'R')
        [phi, D_Q] = Phi.get_phi_2_R(ras_tools, MU, A, ALPHA, sigma, R, N, c_vartheta);
        R2 = R^2 - phi;
        R = sqrt(R2); 
        D_Q_g(i) = D_Q;
    elseif strcmp(ITER_METHOD, 'TEST')
        phi = Phi.TEST_FCT_get_phi_R(ras_tools, MU, A, ALPHA, sigma, R, N, c_vartheta);
        R2 = R^2 - phi;
        R = sqrt(R2); 
    elseif strcmp(ITER_METHOD, 'Y')
        phi = Phi.get_phi_2(ras_tools, MU, LAM, A, ALPHA, sigma, y, e_11, e_20);
        %phi = Phi.get_phi_2_lgPop(ras_tools, MU, LAM, A, ALPHA, sigma, y);
        y2 = y.^2 - phi;
        y = sqrt(y2);  
    elseif strcmp(ITER_METHOD, 'SIM')
        %phi = Phi.phi_2_i_overNcomp(ras_tools, A, ALPHA, MU, LAM, sigma, y, {'PARFOR'});
        phi = Phi.phi_2_experiment(1e4, FIT, MU, LAM, y, sigma);
        y2 = y.^2 - phi;
        y = sqrt(y2);        
    elseif strcmp(ITER_METHOD, 'PHI1')
        %% Phi1 via largePop
        bool_exp = 1; label_di = '$\exp[..]d_i$';
        [phi, ki, di] = Phi.get_phi_1(ras_tools, MU, LAM, A, ALPHA, sigma, y, bool_exp);
        %% Phi1 via C and D_Q
        %bool_di = 0; label_di = '$d_i$';
        %[phi, ki, di] = Phi.get_phi_1_viaC_DQ(ras_tools, C_MU_LAM, A, ALPHA, sigma, y, bool_di);
        %% Phi1 via C and D_i
        %bool_di = 1; label_di = '$d_i$';
        %[phi, ki, di] = Phi.get_phi_1_viaC_Di(ras_tools, C_MU_LAM, A, ALPHA, sigma, y, bool_di);
        %% Iterate
        k_g(:,i) = ki;
        d_g(:,i) = di;
        y = y - phi;       
    else
        error('ITER_METHOD not correctly defined.')
    end
    
    phi_g(:,i) = phi;
    
    if F_g(i) < F_STOP || R_g(i) < R_STOP
        fprintf('F_STOP/R_STOP \n');
        break
    end
    if any(imag(y) ~= 0)
        error('IMAG BREAK');
    end
end

%title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha$=$2\pi$', ', $A$=',num2str(A), ', $N$=',num2str(N),', $\sigma^*$=',num2str(SIGMA_NORM,'%.2f'), ', rng(',num2str(SEED), ')',', IT=',ITER_METHOD,''];

%% Plot iterated dynamics
c = linspecer(4);
% fig = figure; 
    hold on;
    gen = 0:N_G; gen(isnan(R_g))=nan;
    %title(title_str); legend('location','best'); xlabel('$g$'); ylabel('Iter. Dynamics');
    plot(gen, R_g, '-.', 'color', 'r', 'DisplayName', ['Iter: ',ITER_METHOD]);
    %plot(gen, F_g, '-.', 'color', c(2,:),'DisplayName', '$F(g)$');
    %plot(gen, sigma_g, '-.', 'color', c(3,:),'DisplayName', '$\sigma(g)$');
    %plot(gen, sigma_g*N./R_g, '-.', 'color', c(4,:),'DisplayName', '$\sigma^*(g)$');
    set(gca, 'YScale', 'log');

end
    