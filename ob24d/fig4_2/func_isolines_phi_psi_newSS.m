function[s_phi,r_phi] = func_isolines_phi_psi_newSS(MU,LAM,N,A,ALPHA,TAU,SHOW_PLOT)


%% Population
% MU = 100;
% LAM = 200;
% N = 100;
% A = 1;
% ALPHA = 2*pi;
% TAU = 1/sqrt(2*N);

e_10 = e_vartheta_a_b(MU/LAM, 1, 0);
e_11 = e_vartheta_a_b(MU/LAM, 1, 1);

sigma_eps = sqrt(N/2)*A;
R_inf = sqrt(sigma_eps*N/(4*e_10*MU));

%% Plot
step = 1;
R_start = log10(R_inf/5);
R_end = log10(R_inf*2);
NUM_VAL = 501;

%% Init values
R_LIST = 10.^linspace(R_start,R_end,NUM_VAL);
% [x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});
x_2ndZero = sqrt( sqrt(N*(8*e_10^2*MU^2 + N)) - N);
SIGMA_NORM_LIST = linspace(0.01, ceil(x_2ndZero), NUM_VAL); %1.2*x_2ndZero

%% Get phi_R
[mesh_R, mesh_sigmaNorm] = meshgrid(R_LIST, SIGMA_NORM_LIST);
mesh_sigma = mesh_sigmaNorm.*mesh_R/N;
mesh_phi_R = nan*mesh_sigma;
mesh_psi_R = nan*mesh_sigma;

%% PSI=0
p = @(sign,R) 1/sign^2+1/(2*N)+N^3*A^2/(8*sign^4*R^4)-4*e_10^2;

for s=1:length(mesh_sigma)
    for r=1:length(mesh_R)
        sigma = mesh_sigma(s,r);
        R = mesh_R(s,r);
        mesh_phi_R(s,r) = phi_R(MU, LAM, N, A, ALPHA, sigma, R);

        %% CHOoSE PSI
        mesh_psi_R(s,r) = psi_R(A, ALPHA, TAU, sigma, R, N, e_10, e_11);
        %mesh_psi_R(s,r) = psi_R_vanishExp(A, ALPHA, TAU, sigma, R, N, e_10, e_11);
        %mesh_psi_R(s,r) = p(sigma/R*N, R);
            
    end
end

mesh_phi_R = mesh_phi_R./(mesh_R.^2)*N/2;

max_phi = ceil(max(mesh_phi_R,[],'all'));
min_phi = floor(min(mesh_phi_R,[],'all'));

levels = min_phi:step:max_phi;

%% Plot ISO-lines
if SHOW_PLOT==1
    fig1 = figure; hold on;
        colormap('parula'); 
        cbar = colorbar;
        ylabel(cbar, '$\varphi_R^{\mathrm{II},*}$');
        [cfig,h] = contourf(mesh_sigmaNorm, mesh_R, mesh_phi_R, 'LevelList',levels);
    title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha=2\pi', ', A=',num2str(A), ', N=',num2str(N), '$'];
    title(title_str); 
    xlabel('$\sigma^*$'); ylabel('$R$');
    xlim([0,SIGMA_NORM_LIST(end)]);
    ylim([0,R_LIST(end)]);
    set(gca, 'YScale', 'log');
    myfigstyle(gcf, 16, 10, 9, 9);

    %% Plot line of no progress
    line = 'w';
    contour(mesh_sigmaNorm, mesh_R, mesh_phi_R,[0 0],line,'linew',3)
    ylim([R_LIST(1),R_LIST(end)])
    yline(R_inf, 'k--')
end

%% PLOT SETTINGS
% xticks([0:10:50]);
% yticks([1e-1,1,1e1,1e2]);
% set(cbar,'YTick',[-8:2:6])

%% SAR ISO
% fig2 = figure; hold on;
%     colormap('cool'); 
%     cbar = colorbar;
%     ylabel(cbar, '$\psi$');
%     [cfig,h] = contourf(mesh_sigmaNorm, mesh_R, mesh_psi_R, 9);
% % title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha=2\pi', ', A=',num2str(A), ', N=',num2str(N), '$'];
% % title(title_str); 
% xlabel('$\sigma^*$'); ylabel('$R$');
% xlim([0,SIGMA_NORM_LIST(end)]);
% ylim([0,R_LIST(end)]);
% set(gca, 'YScale', 'log');
% myfigstyle(gcf, 16, 10, 9, 9);
% 
% %% Plot line of no progress
% contour(mesh_sigmaNorm, mesh_R, mesh_psi_R,[0 0],'w','linew',3)
% ylim([R_LIST(1),R_LIST(end)])
% % xticks([0:10:50]);
% % yticks([1e-1,1,1e1,1e2]);
% % Plot constant sigma
% sigma0 = (A/2)^(1/2)*1/(8*e_10^2*N-1)^(1/4);
% SIGMA_NORM = sigma0*N./R_LIST;
% plot(SIGMA_NORM, R_LIST, 'k--','linew', 2);
% 
% %% BACK TO PROGRESS RATE
% figure(fig1); hold on;
% contour(mesh_sigmaNorm, mesh_R, mesh_psi_R,[0 0],'m--','linew',3)
% plot(SIGMA_NORM, R_LIST, 'k--','linew', 2);
% 
% % R0
% sigma_eps = sqrt(N/2)*A;
% R_inf = sqrt(sigma_eps*N/(4*e_10*MU));
% R0_v1 = (N*A/2/sqrt(8*e_10^2*N-1)/(sqrt(1+8*e_10^2*MU^2/N)-1))^0.5;
% R0_v2 = N*A/(32*e_10^2*MU);
% R_ratio = R_inf/R0_v1;
% 
% c = e_10;
% m = MU;
% s_0 = sigma0;
% x1 = sqrt(sqrt(128*c^2*m^2*N*(A^2/s_0^4 + 4) + 64*N^2)/(2*(A^2/s_0^4 + 4)) - (4*N)/(A^2/s_0^4 + 4));
% R1 = sigma0*N/x1;
% plot(x1, R1, 'ko');

%% GET PHI-DIP
filter=prod(mesh_phi_R>0,2);
ids = find(filter==1);
id_sig = ids(end);
[phi_min,id_r] = min(mesh_phi_R(id_sig,:));
s_phi = mesh_sigmaNorm(id_sig,1);
r_phi = mesh_R(1,id_r);
% if SHOW_PLOT==1
%     plot(s_phi,r_phi, 'rx', 'MarkerSize',10);
% end

%% GET PSI-DIP
% filter=prod(mesh_psi_R<0,2);
% ids = find(filter==1);
% id_sig = ids(1);
% [psi_max,id_r] = max(mesh_psi_R(id_sig,:));
% s_psi = mesh_sigmaNorm(id_sig,1);
% r_psi = mesh_R(1,id_r);
% plot(s_psi,r_psi, 'bx', 'MarkerSize',10);

