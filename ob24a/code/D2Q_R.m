function [D2Q, dD2Q_ds] = D2Q_R(A, ALPHA, sigma, R, N)
    f = @(s,r,a,freq,n) 2*n*s.^4 + 4*r.^2.*s.^2 + ...
            n*a^2/2 * (1-exp(-(freq*s).^2)) .* (1-exp(-freq^2*(s.^2 + 2*(r/sqrt(n)).^2))) + ...
            2*n*a*freq^2*s.^2 .* exp(-0.5*freq^2*(s.^2+(r/sqrt(n)).^2)) .* (s.^2 + 2*(r/sqrt(n)).^2);
    df = @(s,r,a,freq,n) 8*n*s^3 + 8*r^2*s + 4*a*freq^2*n*s^3*exp(-(freq^2*(r^2/n + s^2))/2) - a^2*freq^2*n*s*exp(-freq^2*s^2)*(exp(-freq^2*((2*r^2)/n + s^2)) - 1) - a^2*freq^2*n*s*exp(-freq^2*((2*r^2)/n + s^2))*(exp(-freq^2*s^2) - 1) + 4*a*freq^2*n*s*exp(-(freq^2*(r^2/n + s^2))/2)*((2*r^2)/n + s^2) - 2*a*freq^4*n*s^3*exp(-(freq^2*(r^2/n + s^2))/2)*((2*r^2)/n + s^2);

    D2Q = double(f(sigma,R,A,ALPHA,N));
    dD2Q_ds = double(df(sigma,R,A,ALPHA,N));
end