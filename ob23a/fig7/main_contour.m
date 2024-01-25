clear all

sph = Sphere;

%% Population
MU = 100;
LAM = 200;
N = 100;
A = 3;
ALPHA = 2*pi;

%% Plot
step = 1;
R_start = -1;
R_end = 2;
NUM_VAL = 201;

%% Coeff
c_vartheta = e_vartheta_a_b(MU/LAM, 1, 0);
e_11 = e_vartheta_a_b(MU/LAM, 1, 1);

%% Init values
R_LIST = 10.^linspace(R_start,R_end,NUM_VAL);
[x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});
SIGMA_NORM_LIST = linspace(0.01, ceil(1.1*x_2ndZero), NUM_VAL); %1.2*x_2ndZero

%% Get phi_R
[mesh_R, mesh_sigmaNorm] = meshgrid(R_LIST, SIGMA_NORM_LIST);
mesh_sigma = mesh_sigmaNorm.*mesh_R/N;

[mesh_phiR, mesh_D_Q_R] = phi_R(MU, LAM, N, A, ALPHA, mesh_sigma, mesh_R);

mesh_phiR = mesh_phiR./(mesh_R.^2)*N/2;

max_phi = ceil(max(mesh_phiR,[],'all'));
min_phi = floor(min(mesh_phiR,[],'all'));

levels = min_phi:step:max_phi;

%% Plot ISO-lines
figure; hold on;
    colormap('parula'); 
    cbar = colorbar;
    ylabel(cbar, '$\varphi_R^{\mathrm{II},*}$');
    [cfig,h] = contourf(mesh_sigmaNorm, mesh_R, mesh_phiR, 'LevelList',levels);
% title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha=2\pi', ', A=',num2str(A), ', N=',num2str(N), '$'];
% title(title_str); 
xlabel('$\sigma^*$'); ylabel('$R$');

set(gca, 'YScale', 'log');
xlim([0,SIGMA_NORM_LIST(end)]);
ylim([R_LIST(1),R_LIST(end)])

myfigstyle(gcf, 20, 12, 9, 9);

%% Plot line of no progress
line = 'w';
contour(mesh_sigmaNorm, mesh_R, mesh_phiR,[0 0],line,'linew',3)
ylim([0,R_LIST(end)])

%% Plot Sphere phi*=0
% x_zero_analytic = sph.phiZero_medium_analytic(c_vartheta, MU, LAM, N, {'exact'});
% xline(x_opt, 'k:','linew',2);

%% R_phi0 curve
R_zero = @(sign) (N^4*A^2*ones(size(sign))/4./(8*N*c_vartheta^2*MU^2 -2*N*sign.^2 - sign.^4)).^(1/4);

%% Plot approx. zero line using sigma^*(R)
sign_phi_zero = sph.phiZero_medium_analytic(c_vartheta, MU, LAM, N, 'exact')-1e-6;
sign_list = linspace(SIGMA_NORM_LIST(1),sign_phi_zero, 10001);
plot(sign_list, R_zero(sign_list), 'k--','linew', 2);

%% Plot transition relation R(sigma*) from D2_Q transition
sign_list = linspace(SIGMA_NORM_LIST(1),SIGMA_NORM_LIST(end), 10001);
R_trans_DQ = @(SIGMA_NORM, N, ALPHA, A, delta) sqrt(delta*2*N)/ALPHA./sqrt(1+SIGMA_NORM.^2/N);

delta_v1 = 2;
delta_v2 = 5;
sign_list = linspace(SIGMA_NORM_LIST(1),SIGMA_NORM_LIST(end), 10001);

plot(sign_list, R_trans_DQ(sign_list, N, ALPHA, A, delta_v1), 'm-','linew', 2);
plot(sign_list, R_trans_DQ(sign_list, N, ALPHA, A, delta_v2), 'y-','linew', 2);

%% Intersection
c1 = 8*MU^2*c_vartheta^2/N;
c2 = @(delta) ALPHA^4*A^2/16/delta^2;
% v1
sign_int_v1 = sqrt( N*sqrt(1+c1) / sqrt(1+c2(delta_v1)) -N );
plot(sign_int_v1, R_zero(sign_int_v1), 'kx','linew',2);
% v2
sign_int_v2 = sqrt( N*sqrt(1+c1) / sqrt(1+c2(delta_v2)) -N );
plot(sign_int_v2, R_zero(sign_int_v2), 'kx','linew',2);


%% Plot constant sigma
% sigma0 = get_sigma_thresh(ALPHA, A);
% SIGMA_NORM = sigma0*N./R_LIST;
% plot(SIGMA_NORM, R_LIST, 'k--','linew', 2);

% PLOT SETTINGS
myfigstyle(gcf, 8, 5, 9, 9);
xticks([0:10:50]);
yticks([1e-1,1,1e1,1e2]);
set(cbar,'YTick',[-12:4:8])
