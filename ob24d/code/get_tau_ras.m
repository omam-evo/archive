function [TAU] = get_tau_ras(N, ALPHA, A)
    TAU = sqrt(1/N.*(1-8./(ALPHA.^2.*A).*lambertw(0, ALPHA.^3.*A.^(3/2)/8)));
end

