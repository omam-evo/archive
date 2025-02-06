strfitnessfct = 'fsphere'; %'frastrigin10';
mu = 1000;
lam = mu*2;
N = 50;

% R_0 = 1e5;
% v = randn(N,1);
% Y_0 = v/norm(v)*R_0;
Y_0 = 1*ones(N,1);

xmean = Y_0;
stopfitness = 1e-12;
stopeval = 1e6;
sigma = 10;
SET_CMA_ES = 1;
SET_INTERMEDIATE = 1;
SET_MATRIX_NORM = 0;

rng(1);
[xmin] = purecmaes(strfitnessfct, mu, lam, N, xmean, sigma, stopfitness, stopeval, SET_CMA_ES, SET_INTERMEDIATE, SET_MATRIX_NORM);

% figure;
%     plot(vecnorm(p_t,2,1))
%     set(gca, 'YScale', 'log');
%     ylim([0.1,10])