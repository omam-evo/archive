%function [] = isolines_phi_y(MU, LAM, N, A, ALPHA, sig_lim, y_lim, NUM_VAL, sigma_mode, phi_mode)

MU=100;
LAM=200;
N=100;
A=10;
ALPHA=2*pi;
sig_lim=[0,1];
y_lim=[0,5.5];
NUM_VAL=101;
sigma_mode='reg';
phi_mode='gain';

sph = Sphere;

%% CONFIG
% MU = 200;
% LAM = 600;
% N = 100;
% A = 10;
% ALPHA = 2*pi;
% NUM_VAL = 101;

%% Min Max
sig_min = max(sig_lim(1),1e-6);
sig_max = sig_lim(2);
y_min = max(y_lim(1),1e-6);
y_max = y_lim(2);

R_min = y_min*sqrt(N);
R_max = y_max*sqrt(N);

%% CONST SIGMA LINE
% sigma_mode = 'norm';

%% COEFF 
C_MU_LAM = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
c_vartheta = e_vartheta_a_b(MU/LAM, 1, 0);
e_11 = e_vartheta_a_b(MU/LAM, 1, 1);
title_str = ['(',num2str(MU),'/',num2str(MU),', ',num2str(LAM),')-ES', ', $\alpha=2\pi', ', A=',num2str(A), ', N=',num2str(N), '$'];

Y_LIST = linspace(y_min,y_max,NUM_VAL);

%% SIGMA (NORM)
[x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU, LAM, 1, 0), MU, LAM, N, {'full'});
if strcmpi(sigma_mode, 'norm')  
    SIGMA_NORM_LIST = linspace(sig_min, sig_max, NUM_VAL); %1.2*x_2ndZero
    [mesh_Y, mesh_sigmaNorm] = meshgrid(Y_LIST, SIGMA_NORM_LIST);
    mesh_R = mesh_Y*sqrt(N);
    mesh_sigma = mesh_sigmaNorm.*mesh_R/N;   
elseif strcmpi(sigma_mode, 'reg')
    SIGMA_LIST = linspace(sig_min, sig_max, NUM_VAL);
    SIGMA_LIST(1) = 1e-6;
    [mesh_Y, mesh_sigma] = meshgrid(Y_LIST, SIGMA_LIST);
    mesh_R = mesh_Y*sqrt(N);
end

%% NOT NORMALIZED
DQ = @(sigma,R) sqrt(2*N*sigma.^4 + 4*R.^2.*sigma.^2 + ...
            N*A^2/2 * (1-exp(-(ALPHA*sigma).^2)) .* (1-exp(-ALPHA^2*(sigma.^2 + 2*(R/sqrt(N)).^2))) + ...
            2*N*A*ALPHA^2*sigma.^2 .* exp(-0.5*ALPHA^2*(sigma.^2+(R/sqrt(N)).^2)) .* (sigma.^2 + 2*(R/sqrt(N)).^2));
GT = @(sigma, y) 2*y + exp(-0.5*(ALPHA*sigma).^2).*ALPHA*A.*sin(ALPHA*y);

prefac = c_vartheta*mesh_sigma.^2./DQ(mesh_sigma, mesh_R);

%% GET PHI/GAIN
mesh_phi1 = prefac.*GT(mesh_sigma,mesh_Y);
mesh_phi2 = prefac.*(2*mesh_Y.*GT(mesh_sigma,mesh_Y)) - mesh_sigma.^2/MU;
mesh_gain = 2*mesh_Y.*GT(mesh_sigma,mesh_Y);

if strcmpi(phi_mode, 'phi1') 
    mesh_plot = mesh_phi1;
    contour_label_str ='$\varphi_i$';
elseif strcmpi(phi_mode, 'phi2') 
    mesh_plot = mesh_phi2;
    contour_label_str = '$\varphi_i^{\mathrm{II}}$';
elseif strcmpi(phi_mode, 'gain') 
    mesh_plot = mesh_gain;
    
    contour_label_str = 'G($y_i,\sigma)$';
end


%% Plot ISO-lines
figure; hold on;
    %colormap('jet'); 
    cbar = colorbar;
    
    
if strcmpi(sigma_mode, 'norm') 
    %step = 1;
    %mesh_phiR = mesh_phiR./(mesh_R.^2)*N/2;
    %max_phi = ceil(max(mesh_phiR,[],'all'));
    %min_phi = floor(min(mesh_phiR,[],'all'));
    %levels = min_phi:step:max_phi;
    
    [cfig,h] = contourf(mesh_sigmaNorm, mesh_Y, mesh_plot, 25);
    contour(mesh_sigmaNorm, mesh_Y, mesh_plot,[0 0],'w','linew',2)
    %title(title_str); 
    xlabel('$\sigma^*$'); ylabel('$y_i$');
    ylabel(cbar, contour_label_str);
    xlim([sig_min,sig_max]);
    ylim([y_min,y_max]);
    
elseif strcmpi(sigma_mode, 'reg')
    %phi_max = sph.phi(x_opt, C_MU_LAM, MU, LAM, N, 'full');
    levels = -2:0.05:2;
    [cfig,h] = contourf(mesh_sigma, mesh_Y, mesh_plot, 25); %, 'LevelList',levels
    contour(mesh_sigma, mesh_Y, mesh_plot,[0 0],'w','linew',2)
    %title(title_str); 
    xlabel('$\sigma$'); ylabel('$y_i$');
    ylabel(cbar, contour_label_str);
    xlim([sig_min,sig_max]);
    ylim([y_min,y_max]);
end
% end
xline(get_sigma_thresh(ALPHA,A), 'linestyle', ':', 'linewidth', 2, 'alpha', 1)

xlim([sig_lim(1),sig_lim(2)])
ylim([y_lim(1),y_lim(2)])
xticks([0:0.2:1])
yticks([0:1:5])
myfigstyle(gcf, 8, 6, 9, 9);

