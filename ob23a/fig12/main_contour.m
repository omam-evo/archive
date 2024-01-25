clear all

sph = Sphere;

%% Population
MU = 400;
LAM = 800;
N = 100;
A = 10;
ALPHA = 2*pi;

%% Plot
step = 1;
R_start = -1;
R_end = 2;
NUM_VAL = 101;

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

%% Plot constant sigma
sigma0 = get_sigma_thresh(ALPHA, A);
SIGMA_NORM = sigma0*N./R_LIST;
plot(SIGMA_NORM, R_LIST, 'r--','linew', 2);

%% PLOT SETTINGS
myfigstyle(gcf, 8, 5, 9, 9);
xticks([0:20:100]);
yticks([1e-1,1,1e1,1e2]);
set(cbar,'YTick',[-10:5:10])
