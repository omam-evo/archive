clear

ras = RasTools;
A = 1; ALPHA = 2*pi;
N = 100;
TRIALS = 100;
%SigN = 0.3;
SIGMA_NORM = 30;%SigN*N;

R_LIST = 10.^linspace(-3,2,1000)';
D2_Q_exactR = nan*R_LIST;
D2_Q_exactR_2 = nan*R_LIST;
D2_Q_fctR_EV = nan*R_LIST;
D2_Q_fctR = nan*R_LIST;
D2_Q_mean = nan*R_LIST;
D2_Q_stderr = nan*R_LIST;
NormY_mean = nan*R_LIST;
NormY_stderr = nan*R_LIST;

rng(1);
for i=1:length(R_LIST)
    R = R_LIST(i);

    sigma_transf = SIGMA_NORM*R/N;
    
    %% Choosing y_i random with constant R (N-1 dof)
    % random number rescaled to be on hypersphere => jagged curve for small N expected
    v = randn(N-1, 1); 
    Y = v/norm(v)*R;
    D2_Q_exactR(i) = sum(ras.variance_q_gain_simplified(A, ALPHA, Y, sigma_transf));
    
end

figure; 
hold on; legend; title(['$\sigma^*=',num2str(SIGMA_NORM),', N=',num2str(N), '$']);
xlabel('$R$'); ylabel('$D_Q^2(R)$');

plot(R_LIST, D2_Q_exactR, '-', 'DisplayName', 'Exact 1');

set(findall(gcf, '-property', 'XScale'), 'XScale', 'log');
set(findall(gcf, '-property', 'YScale'), 'YScale', 'log');
myfigstyle(gcf, 16, 10, 10, 10);

% subplot(1,2,2); hold on;
% plot(R_LIST, NormY_mean);
% %errorbar(R_LIST, mean_norm_Y, stderr_norm_Y);
% set(findall(gcf, '-property', 'XScale'), 'XScale', 'log');
% set(findall(gcf, '-property', 'YScale'), 'YScale', 'log');
    