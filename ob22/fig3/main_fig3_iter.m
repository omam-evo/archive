clear

%% CONFIG 1a
MU = 150;
LAM = 300;
SET_di = 0; % SET_di = 0: iteration using only k_i-term. SET_di = 1: iteration using f'_i = k_i+d_i
R_0 = 100;

%% CONFIG 1c
% MU = 150;
% LAM = 300;
% SET_di = 1; % SET_di = 1: iteration using f'_i = k_i+d_i
% R_0 = 100;

%% CONFIG 1c
% MU = 150;
% LAM = 300;
% SET_di = 1; % SET_di = 1: iteration using f'_i = k_i+d_i
% R_0 = 0.1; 

%% CONFIG 2a
% MU = 1500;
% LAM = 3000;
% SET_di = 0; 
% R_0 = 100; 

%% CONFIG 2b
% MU = 1500;
% LAM = 3000;
% SET_di = 1; 
% R_0 = 100; 

% %% CONFIG 2c
% MU = 1500;
% LAM = 3000;
% SET_di = 1; 
% R_0 = 0.1; 


%% PARAMS
N = 100;
A = 1;
ALPHA = 2*pi;
SIGMA_NORM = 30;
F_STOP = 1e-6;
N_G = 200;
SEED = 1;
BOOL_PHI_2 = 0;

%% set color
if SET_di==1
    mycolor = 'r--';
elseif SET_di==0
    mycolor = 'b-.';
end

%% Additional plot for R_0=0.1
% Read id from avg. dynamics when R~0.1
% id_read = 75;
% ids = find(~isnan(R_g));
% x = id_read:1:(ids(end)+id_read-1);
% hold on; plot(x, R_g(ids));

%% INIT
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1);
ras_tools = RasTools;
Phi = Phi_v2;
theta = MU/LAM;
try C_MU_LAM = e_mu_lam_a_b(MU, LAM, 1, 0); catch, C_MU_LAM = 1/sqrt(2*pi)/theta * exp(-1/2*(norminv(theta))^2); end 
[order_y_min, order_y_max] = ras_tools.get_all_extrema(ALPHA, A, N);

%% Position
rng(SEED); v = randn(N,1); y = v/norm(v)*R_0;
% y = -1.75*ones(N, 1);

%% Data
y_g = nan*zeros(N, N_G+1);
k_g = nan*zeros(N, N_G+1);
d_g = nan*zeros(N, N_G+1);
phi_g = nan*zeros(N, N_G+1);
R_g = nan*zeros(N_G+1, 1);
F_g = R_g;
sigma_g = R_g;
title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $A$=',num2str(A),', $N$=',num2str(N), ', rng(',num2str(SEED), ')'];

for i=1:N_G+1
    R = norm(y);
    sigma = SIGMA_NORM*R/N;
    
    R_g(i) = R;
    F_g(i) = FIT(y);
    sigma_g(i) = sigma;
    y_g(:,i) = y;
    
    %% Try numerical exact
%     y_i = y(1); y_inj = y(2:end);
%     phi_test = Phi(MU, LAM, ras_tools, A, ALPHA, y, y_i, y_inj, sigma, 'nonlinear');
%     phi_i = phi_test.get_phi_genericPQ_num();
%     phi_g = phi_i*ones(N,1);
    
    fprintf('\t Iter %i \t R(g)=%d \t F(g)=%d \t sigma(g)=%d \n', i, R_g(i), F_g(i), sigma_g(i));
    if BOOL_PHI_2 == 0
        %% Phi1 via largePop
%         bool_exp = 1; label_di = '$\exp[..]d_i$';
%         [phi, ki, di] = Phi.get_phi_1(ras_tools, MU, LAM, A, ALPHA, sigma, y, bool_exp);
        %% Phi1 via C and D_Q
%         bool_di = 0; label_di = '$d_i$';
%         [phi, ki, di] = Phi.get_phi_1_viaC_DQ(ras_tools, C_MU_LAM, A, ALPHA, sigma, y, bool_di);
        %% Phi1 via C and D_i
        bool_di = SET_di; label_di = '$d_i$';
        [phi, ki, di] = Phi.get_phi_1_viaC_Di(ras_tools, C_MU_LAM, A, ALPHA, sigma, y, bool_di);
        %% Iterate
        k_g(:,i) = ki;
        d_g(:,i) = di;
        y = y - phi;
    else
        bool_phi2_exp = 1; label_di = '$\varphi^{II}$';
        phi = Phi.get_phi_2(ras_tools, MU, LAM, A, ALPHA, sigma, y, bool_phi2_exp);
        y_sq = y.^2 - phi;
        y = sqrt(y_sq);
    end
    phi_g(:,i) = phi;

    if F_g(i) < F_STOP
        fprintf('F_STOP \n');
        break
    end
    if any(imag(phi) ~= 0)
        fprintf('IMAG BREAK \n');
        break
    end
end

fig = figure; hold on;
    gen = 0:N_G; gen(isnan(R_g))=nan;
    title(title_str); legend; xlabel('$g$'); ylabel('Iter. dynamics');
    plot(gen, R_g, mycolor, 'DisplayName', '$R(g)$');
    %plot(gen, F_g, 'g-', 'DisplayName', '$F(g)$');
    %plot(gen, sigma_g, 'b-', 'DisplayName', '$\sigma(g)$');
    %plot(gen, sigma_g*N./R_g, 'c-', 'DisplayName', '$\sigma^*(g)$');
    set(gca, 'YScale', 'log');

myfigstyle(fig, 30, 15, 12, 12);
