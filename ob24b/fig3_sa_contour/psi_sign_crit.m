%% PARAMS
MU = 100;
LAM = 200; 
N = 20; 
A = 1;
R_inf = sqrt(sqrt(N/2)*A*N/(4*e_10*MU));

TAU = 1/sqrt(2*N);

e_10 = e_vartheta_a_b(MU/LAM, 1, 0);
e_11 = e_vartheta_a_b(MU/LAM, 1, 1);


f1 = @(x) 2*e_10*x*(1+0.5*ALPHA^2*A*exp(-(ALPHA*x)^2*A/(sqrt(128*N)*e_10*MU))) ...
    -sqrt(1 + x^2/(2*N) + 4*e_10^2*MU^2/x^2);

x0 = 0.1; %sqrt( sqrt(N*(8*e_10^2*MU^2 + N)) - N); %sqrt(MU);
fzero(f1,x0)

f_psi = @(x) psi_R(A, ALPHA, TAU, x*R_inf/N, R_inf, N, e_10, e_11);
fzero(f_psi,x0)