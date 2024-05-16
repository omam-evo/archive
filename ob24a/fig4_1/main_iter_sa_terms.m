clear all
addpath('../n0_generic')

%% CONFIG
NEGLECT_TERMS = 0;

MU = 100;
LAM = 200;
A = 1;
N = 100;
N_G = 10000;
TAU = 1/sqrt(8*N); %GAMMA = sqrt(1-1/8);
ITER_METHOD = 'R'; % 'SIM' using simulations, 'R' via phi^(II)(R),'Y' via phi^(II)

ALPHA = 2*pi;
F_STOP = 1e-2; 
R_STOP = nan;  
SIGMA_STOP = 1e-4;


%% Coeff.
e_10 = e_vartheta_a_b(MU/LAM, 1, 0);
e_11 = e_vartheta_a_b(MU/LAM, 1, 1); % Adjust to include terms for phi^(II)(y_i)
e_20 = e_vartheta_a_b(MU/LAM, 2, 0); % Adjust to include terms for phi^(II)(y_i)
e_11_setPhi2 = 0;
e_20_setPhi2 = 0;

%% INIT
FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1);
FIT_R = @(R) R^2 + N*A*(1 - exp(-0.5*(ALPHA*R/sqrt(N))^2));
ras_tools = RasTools;
sph = Sphere;
phi_class = Phi_v2;
psi_class = Psi;
[order_y_min, order_y_max] = ras_tools.get_all_extrema(ALPHA, A, N);

%% Position and Mutation
% a = 20;   % R > 7*sqrt(N) for A=10 necessary
% R_0 = a*sqrt(N);
% % rng(SEED); 
% v = randn(N,1); 
% y = v/norm(v)*R_0;
y = 50*ones(N,1); 
R_0 = norm(y);
x_2ndZero = sqrt( sqrt(N*(8*e_10^2*MU^2 + N)) - N);
sigma = x_2ndZero * R_0 / N;        

%% Data
y_g = nan*zeros(N, N_G+1);
k_g = nan*zeros(N, N_G+1);
d_g = nan*zeros(N, N_G+1);
R_g = nan*zeros(N_G+1, 1);
F_g = nan*zeros(N_G+1, 1);
sigma_g = nan*zeros(N_G+1, 1);
D_Q_g = nan*zeros(N, N_G+1);
psi_g = nan*zeros(N_G+1, 1);

if strcmp(ITER_METHOD, 'R') || strcmp(ITER_METHOD, 'TEST')
    R = R_0;
    phi_g = nan*zeros(1, N_G+1);
else
    phi_g = nan*zeros(N, N_G+1);
end

for i=1:N_G+1
    
    if strcmp(ITER_METHOD, 'R')   
        F_g(i) = FIT_R(R);
    else %'SIM', 'Y'
        R = norm(y);
        F_g(i) = FIT(y); 
        y_g(:,i) = y;
    end
    
    %% Save pos. and mut. now because overwritten during iteration
    R_g(i) = R;
    sigma_g(i) = sigma;
    
    fprintf('\t Iter %i, R(g)=%d, F(g)=%d, sigma(g)=%d \n', i, R_g(i), F_g(i), sigma_g(i));

    if strcmp(ITER_METHOD, 'R')
        if NEGLECT_TERMS == 0       
            [phi, ~] = phi_class.get_phi_2_R(ras_tools, MU, A, ALPHA, sigma, R, N, e_10);
            [psi,~,~,~] = psi_class.get_psi_v1_R(ras_tools, A, ALPHA, TAU, sigma, R, N, e_10, e_11);
        else
            phi = phi_R_vanishExp(MU, LAM, N, A, ALPHA, sigma, R);
            psi = psi_R_vanishExp(A, ALPHA, TAU, sigma, R, N, e_10, e_11);
        end
        % update position
        R2 = R^2 - phi;
        R = sqrt(R2); 
        % update sigma
        sigma = sigma*(1+psi);
    elseif strcmp(ITER_METHOD, 'SIM')
        [phi,~,psi,~] = phi_class.phi2_psi_experiment(1e3, FIT, MU, LAM, TAU, y, sigma);
        y2 = y.^2 - phi;
        y = sqrt(y2);    
        sigma = sigma*(1+psi);
    else
        error('ITER_METHOD not correctly defined.')
    end
    
    % Save
    phi_g(:,i) = phi;
    psi_g(i) = psi;
    
    if F_g(i) < F_STOP || R_g(i) < R_STOP
        fprintf('F_STOP/R_STOP \n');
        break
    end
    if sigma_g(i) < SIGMA_STOP
        fprintf('SIGMA_STOP \n');
        break
    end
    if any(imag(y) ~= 0)
        error('IMAG BREAK');
    end
end

title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha$=$2\pi$', ', $A$=',num2str(A), ', $N$=',num2str(N),', $\tau$=',num2str(TAU,'%.3f'),', IT=',ITER_METHOD,'']; %, ', rng(',num2str(SEED), ')'

%% Plot iterated dynamics
if strcmpi(ITER_METHOD, 'Y')
    myc = 'b-.';
elseif strcmpi(ITER_METHOD, 'R')
    myc = 'r--';
elseif strcmpi(ITER_METHOD, 'SIM')
    myc = 'g-.';
else
    myc = 'm-.';
end
% fig = figure; 
%     hold on;
%     gen = 0:N_G; gen(isnan(R_g))=nan;
%     title(title_str); legend('location','best'); xlabel('$g$'); ylabel('Iter. Dynamics');
%     plot(gen, R_g, myc, 'DisplayName', ['Iter. ',ITER_METHOD]);
%     plot(gen, sigma_g*N./R_g, 'm:', 'DisplayName', ['$\sigma^*$']);
%     set(gca, 'YScale', 'log');
%     myfigstyle(fig, 16, 10, 10, 10);  
sigma_eps = sqrt(N/2)*A;
R_inf = sqrt(sigma_eps*N/(4*e_10*MU));
sign_crit = get_sign_w(N, A, ALPHA, MU, LAM);

fig = figure; 
    hold on;
    title(title_str); legend('location','best'); xlabel('$g$'); ylabel('Iter. Dynamics');
    plot(sigma_g*N./R_g, R_g, myc, 'DisplayName', ['Iter. ',ITER_METHOD]);
    sigma_eps = sqrt(N/2)*A;
    %yline(sqrt(sigma_eps*N/(4*e_10*MU)), 'k--');
    %plot(sign_crit, R_inf, 'o');
    set(gca, 'YScale', 'log');
    myfigstyle(fig, 16, 10, 10, 10);  
   