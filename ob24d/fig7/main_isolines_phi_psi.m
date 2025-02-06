clear all
addpath('../n0_generic')

sph = Sphere;

%% Population
MU = 100;
LAM = 200;
N = 100;
A = 1;
ALPHA = 2*pi;
TAU = 1/sqrt(8*N);

%% Plot
step = 1;
R_start = -1;
R_end = 2;
NUM_VAL = 201;

%% Coeff
e_10 = e_vartheta_a_b(MU/LAM, 1, 0);
e_11 = e_vartheta_a_b(MU/LAM, 1, 1);

%% Init values
R_LIST = 10.^linspace(R_start,R_end,NUM_VAL);
% [x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});
x_2ndZero = sqrt( sqrt(N*(8*e_10^2*MU^2 + N)) - N);
SIGMA_NORM_LIST = linspace(0.01, ceil(1.1*x_2ndZero), NUM_VAL); %1.2*x_2ndZero

%% Get phi_R
[mesh_R, mesh_sigmaNorm] = meshgrid(R_LIST, SIGMA_NORM_LIST);
mesh_sigma = mesh_sigmaNorm.*mesh_R/N;
mesh_phi_R = nan*mesh_sigma;
mesh_psi_R = nan*mesh_sigma;

for s=1:length(mesh_sigma)
    for r=1:length(mesh_R)
        sigma = mesh_sigma(s,r);
        R = mesh_R(s,r);
        mesh_phi_R(s,r) = phi_R(MU, LAM, N, A, ALPHA, sigma, R);
%         if strcmpi(set_pos_mode, 'R')
            mesh_psi_R(s,r) = psi_R(A, ALPHA, TAU, sigma, R, N, e_10, e_11);
%         elseif strcmpi(set_pos_mode, 'Y')
%             Y = R/sqrt(N)*ones(N,1);
%             mesh_psi_R(s,r) = psi_Y(A, ALPHA, TAU, sigma, Y, e_10, e_11);
%             %mesh_psi_R(s,r) = psi_hybrid(A, ALPHA, TAU, sigma, Y, e_10, e_11);
%         end
    end
end

mesh_phi_R = mesh_phi_R./(mesh_R.^2)*N/2;

max_phi = ceil(max(mesh_phi_R,[],'all'));
min_phi = floor(min(mesh_phi_R,[],'all'));

levels = min_phi:step:max_phi;

%% Plot ISO-lines
figure; hold on;
    colormap('parula'); 
    cbar = colorbar;
    ylabel(cbar, '$\varphi_R^{\mathrm{II},*}$');
    [cfig,h] = contourf(mesh_sigmaNorm, mesh_R, mesh_phi_R, 'LevelList',levels);
% title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha=2\pi', ', A=',num2str(A), ', N=',num2str(N), '$'];
% title(title_str); 
%xlabel('$\sigma^*$'); 
ylabel('$R$');
xlim([0,SIGMA_NORM_LIST(end)]);
ylim([0,R_LIST(end)]);
set(gca, 'YScale', 'log');
myfigstyle(gcf, 8, 3.5, 8, 8);

%% Plot line of no progress
line = 'w';
contour(mesh_sigmaNorm, mesh_R, mesh_phi_R,[0 0],line,'linew',3)
ylim([R_LIST(1),R_LIST(end)])

%% PLOT SETTINGS
xticks([0:10:50]);
yticks([1e-1,1,1e1,1e2]);
set(cbar,'YTick',[-8:2:6])

%% Plot constant sigma
% sigma0 = get_sigma_thresh(ALPHA, A);
% SIGMA_NORM = sigma0*N./R_LIST;
% plot(SIGMA_NORM, R_LIST, 'k--','linew', 2);

%% SAR ISO
figure; hold on;
    colormap('cool'); 
    cbar = colorbar;
    ylabel(cbar, '$\psi$');
    [cfig,h] = contourf(mesh_sigmaNorm, mesh_R, mesh_psi_R, 9);
% title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha=2\pi', ', A=',num2str(A), ', N=',num2str(N), '$'];
% title(title_str); 
xlabel('$\sigma^*$'); ylabel('$R$');
xlim([0,SIGMA_NORM_LIST(end)]);
ylim([0,R_LIST(end)]);
set(gca, 'YScale', 'log');
myfigstyle(gcf, 8, 4, 8, 8);

%% Plot line of no progress
line = 'w';
contour(mesh_sigmaNorm, mesh_R, mesh_psi_R,[0 0],line,'linew',3)
ylim([R_LIST(1),R_LIST(end)])
xticks([0:10:50]);
yticks([1e-1,1,1e1,1e2]);