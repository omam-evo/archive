function [res] = e_integrand(x, mu, lam, a, b, logfac)
    
    phi = normcdf(x); % phi in [0,1]
    if lam-mu-1 ~= 0
        logcdf = (lam-mu-1)*log(phi);  %init
        logcdf(phi==0) = 0; %exp[log(0)] = 0
    else
        logcdf = 0*x; %(lam-mu-1)*log(normphi(x)) = 0
    end
    if mu-a ~= 0    % ToDo what if a > MU
        logmincdf = (mu-a)*log(1-phi);  %init
        logmincdf((1-phi)==0) = 0; %exp[log(0)] = 0; problem with (0)^(-1) is mitigated
    else
        logmincdf = 0*x; %(lam-mu-1)*log(normphi(x)) = 0
    end
    
    res = x.^b .* exp( logfac-(a+1)/2*x.^2 + logcdf + logmincdf );
end