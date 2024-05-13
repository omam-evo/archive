%% Analogous to SMN C.23; SAR on general function with R-dependent EQ and D2Q
function [psi,dDQ,DQ,dEQ] = psi_R_vanishExp(A, ALPHA, TAU, sigma, R, N, e_10, e_11)

    dEQ = @(s,n) 2*n.*s;
    SET_ZERO = 0;
    D2_Q = @(s,r,a,freq,n) 2*n*s.^4 + 4*r.^2.*s.^2 + ...
            n*a^2/2 * (1-SET_ZERO*exp(-(freq*s).^2)) .* (1-SET_ZERO*exp(-freq^2*(s.^2 + 2*(r/sqrt(n)).^2))) + ...
            2*n*a*freq^2*s.^2 .* SET_ZERO * exp(-0.5*freq^2*(s.^2+(r/sqrt(n)).^2)) .* (s.^2 + 2*(r/sqrt(n)).^2);

    psi = TAU^2*(0.5 - e_10*sigma*dEQ(sigma,N)/sqrt(D2_Q(sigma,R,A,ALPHA,N)));  

end