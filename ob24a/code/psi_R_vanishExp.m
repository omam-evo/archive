%% Analogous to SMN C.23; SAR on general function with R-dependent EQ and D2Q
function [psi,dDQ,DQ,dEQ] = psi_R_vanishExp(A, ALPHA, TAU, sigma, R, N, e_10, e_11)

    ZERO = 0;
    dEQ = @(s, r, a, freq, n) 2*n*s + a*freq^2*n*s*exp(-0.5*(freq^2*r^2)/n)*exp(-0.5*(freq^2*s^2));

    psi = TAU^2*(0.5 - e_10*sigma*dEQ(sigma,R,A,ALPHA,N)/sqrt(D2Q_R_terms(A, ALPHA, sigma, R, N)));  

end