function [s, r] = isolines_crit_point(MU, LAM, N, A, ALPHA, NUM_VAL, R_end, SHOW_PLOT, STEP_PLOT, sigma0)

% phi = Phi_v2;
% ras_tools = RasTools;
sph = Sphere;
iter = IterVariants;

%% CONFIG
% MU = 100;
% LAM = 200;
% N = 100;
% A = 1;
% step = 1;
% R_end = 20;

% ALPHA = 2*pi;
C_MU_LAM = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
c_vartheta = e_vartheta_a_b(MU/LAM, 1, 0);
e_11 = e_vartheta_a_b(MU/LAM, 1, 1);
title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha=2\pi', ', A=',num2str(A), ', N=',num2str(N), '$'];

R_start = 0;

% NUM_VAL = 201;
R_LIST = linspace(R_start,R_end,NUM_VAL);
R_LIST(1) = 1e-4; %10.^linspace(-2,2,100);
[x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});
SIGMA_NORM_LIST = linspace(0.01, ceil(1.1*x_2ndZero), NUM_VAL); %1.2*x_2ndZero

[mesh_R, mesh_sigmaNorm] = meshgrid(R_LIST, SIGMA_NORM_LIST);

mesh_sigma = mesh_sigmaNorm.*mesh_R/N;

%% NOT NORMALIZED
[mesh_phiR, mesh_D_Q_R] = iter.get_phi_R(MU, A, ALPHA, mesh_sigma, mesh_R, N, c_vartheta);

mesh_phiR = mesh_phiR./(mesh_R.^2)*N/2;

max_phi = ceil(max(mesh_phiR,[],'all'));
min_phi = floor(min(mesh_phiR,[],'all'));

levels = min_phi:STEP_PLOT:max_phi; %min_phi:step:max_phi; %needed for v4: levels = -14:1:6

%% Plot ISO-lines
if SHOW_PLOT
    figure; hold on;
        colormap('jet'); 
        cbar = colorbar;
        ylabel(cbar, '$\varphi_R^{\mathrm{II},*}$');
        [cfig,h] = contourf(mesh_sigmaNorm, mesh_R, mesh_phiR, 'LevelList',levels);
        title(title_str); xlabel('$\sigma^*$'); ylabel('$R$');
    xlim([0,SIGMA_NORM_LIST(end)])
    ylim([0,R_LIST(end)])
    myfigstyle(gcf, 20, 12, 9, 9);

    %% Plot line of no progress
    line = 'w';
    contour(mesh_sigmaNorm, mesh_R, mesh_phiR,[0 0],line,'linew',3)
    ylim([0,R_LIST(end)]);
end

%% Fid min value
% [phi,R,signorm]=isolines(100, 200, 100, 1, 2*pi, 1, 20);
filter=prod(mesh_phiR>0,2);
ids = find(filter==1);
id_sig = ids(end);
[phi_min,id_r] = min(mesh_phiR(id_sig,:));
s = mesh_sigmaNorm(id_sig,1);
r = mesh_R(1,id_r);
if SHOW_PLOT
    plot(s,r, 'ko');
end

%% Plot Sphere phi*=0
% x_zero_analytic = sph.phiZero_medium_analytic(c_vartheta, MU, LAM, N, {'exact'});
% xline(x_opt, 'k:','linew',2);

%% Plot approx. zero line using sigma^*(R)
% Rstart_sign_zero = (N^3*A^2/(32*c_vartheta^2*MU^2))^(1/4);
% Rend_sign_zero = max(R_LIST); 
% eps = 1e-6;
% R_LIST_sign_zero = linspace(Rstart_sign_zero+eps, Rend_sign_zero, NUM_VAL);
% sign_zero = ( (N^2 - N^4*A^2/4*ones(size(R_LIST_sign_zero))./(R_LIST_sign_zero.^4) + 8*N*c_vartheta^2*MU^2).^(1/2) - N ).^(1/2);
% plot(sign_zero, R_LIST_sign_zero, 'k--','linew', 2);

%% Plot transition relation R(sigma*) from D2_Q transition
% delta = 5;
% R_trans_DQ = @(SIGMA_NORM, N, ALPHA, A) sqrt(delta*2*N)/ALPHA./sqrt(1+SIGMA_NORM.^2/N);
% plot(SIGMA_NORM_LIST, R_trans_DQ(SIGMA_NORM_LIST, N, ALPHA, A), 'm--','linew', 2);

%% Plot requirement sigma*>...
% R < alpha*A*N/2
% sigma_bound = @(N,ALPHA,A,R) sqrt(2)*N/ALPHA*sqrt(log(ALPHA*A*N/2./R))./R;
% sigma_bound_v2 = @(N,ALPHA,A,R) sqrt(2)*N/ALPHA*sqrt(log(ALPHA*A*sqrt(N)/2./R))./R;
% plot(sigma_bound(N,ALPHA,A,R_LIST), R_LIST, 'c--','linew', 2);
% plot(sigma_bound_v2(N,ALPHA,A,R_LIST), R_LIST, 'b--','linew', 2);

%% Tests
% sigma_bound_v3 = @(N,ALPHA,A,R) sqrt(2)*N/ALPHA*sqrt(log(3/8*ALPHA*A*N./(R.^2)))./R;
% sigma_bound_v4 = @(N,ALPHA,A,R) sqrt(2)*N/ALPHA*sqrt(log( ALPHA^4*A*N/12./(R.^2)-ALPHA^2*A/2 ))./R;
% R_bound_v1 = @(N,ALPHA,A,signorm) 1/sqrt(2)*sqrt( sqrt(N^4*(2+ALPHA^2*A)^2/ALPHA^4./(signorm.^4) + 4*ALPHA^2*A*N^3/6./(signorm.^2)) - N^2*(2+ALPHA^2*A)/ALPHA^2./(signorm.^2) );
% plot(sigma_bound_v3(N,ALPHA,A,R_LIST), R_LIST, 'g--','linew', 2);
% plot(sigma_bound_v4(N,ALPHA,A,R_LIST), R_LIST, 'y--','linew', 2);
% plot(SIGMA_NORM_LIST, R_bound_v1(N,ALPHA,A,SIGMA_NORM_LIST) , 'y--','linew', 2);

%% Plot constant sigma
SIGMA_NORM = sigma0*N./R_LIST;
if SHOW_PLOT
    plot(SIGMA_NORM, R_LIST, 'k--','linew', 2);
end

%% Crossing line sigma*const and phi=0
% signorm_sec = sqrt( (sqrt((1+A^2/4/sigma0^4)*8*N*c_vartheta^2*MU^2+N^2)-N)/(1+A^2/4/sigma0^4));
% plot(signorm_sec, sigma0*N/signorm_sec, 'ko', 'MarkerSize',15);

% %% SAR-Isolines
% TAU = 1/sqrt(2*N);
% psi_class = Psi;
% mesh_psi = nan*mesh_sigma;
% for i=1:length(R_LIST)
%     for j=1:length(SIGMA_NORM_LIST)
%         mesh_psi(i,j) = psi_class.get_psi_v1_R(ras_tools, A, ALPHA, TAU, mesh_sigma(i,j), mesh_R(i,j), N, c_vartheta, e_11);
%     end
% end
% 
% % mesh_phiR = real(mesh_phiR);warning('REAL')
% C = mesh_phiR; 
% C(C<=0) = -inf;
% 
% max_phi = ceil(max(mesh_phiR,[],'all'));
% min_phi = floor(min(mesh_phiR,[],'all'));
% 
% levels = min_phi:step:max_phi; %min_phi:step:max_phi; %needed for v4: levels = -14:1:6
% 
% %% Plot SAR-ISO-lines
% max_psi = max(mesh_psi,[],'all');
% min_psi = min(mesh_psi,[],'all');
% %% LEVELS: add =>   , 'LevelList',levels
% steps = 15;
% step = (max_psi-min_psi)/steps;
% levels_plus = [0:step:max_psi, max_psi]; 
% levels_minus = -fliplr([step:step:-min_psi, -min_psi]); 
% levels = [levels_minus,levels_plus];
% figure; hold on;
%     colormap('jet'); 
%     cbar = colorbar;
%     ylabel(cbar, '$\psi_R$');
%     [cfig,h] = contourf(mesh_sigmaNorm, mesh_R, mesh_psi, 'LevelList',levels);
%     title([title_str, ', $\tau$=',num2str(TAU)]); xlabel('$\sigma^*$'); ylabel('$R$');
% xlim([0,SIGMA_NORM_LIST(end)])
% ylim([0,R_LIST(end)])
% myfigstyle(gcf, 12, 7, 9, 9);
% 
% %% Plot line of SAR=0
% hold on
% contour(mesh_sigmaNorm, mesh_R, mesh_psi, [0 0], 'k-', 'linew',3)
% ylim([0,R_LIST(end)])
end