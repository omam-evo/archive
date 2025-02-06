function [EQ, dEQ_ds] = EQ_R(A, ALPHA, sigma, R, N)
    f = @(s, r, a, freq, n) n*s.^2 + a*n*exp(-0.5*freq^2*r^2/N).*(1-exp(-0.5*freq^2*s.^2));
    df = @(s, r, a, freq, n) 2*n*s + a*freq^2*n*s*exp(-0.5*(freq^2*r^2)/n)*exp(-0.5*(freq^2*s^2));

    EQ = double(f(sigma,R,A,ALPHA,N));
    dEQ_ds = double(df(sigma,R,A,ALPHA,N));
end