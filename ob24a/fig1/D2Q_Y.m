function [D2Q, dD2Q_ds] = D2Q_Y(A, ALPHA, sigma, Y)
    fi = @(s,y,A,ALPHA) 4*Y.^2.*s.^2 + 2*s.^4 + ...
        A^2/2.*(1 - exp(-(ALPHA*s).^2)).*(1-cos(2*ALPHA*Y).*exp(-(ALPHA*s).^2)) + ...
        2*A*ALPHA*sigma.^2.*exp(-0.5*(ALPHA*s).^2).*(ALPHA*s.^2.*cos(ALPHA*Y) + 2*Y.*sin(ALPHA*Y));
    dfi_ds = @(s,y,A,ALPHA) s.*y.^2.*8.0+s.^3.*8.0-A.*ALPHA.^3.*s.^3.*exp(ALPHA.^2.*s.^2.*(-1.0./2.0)).*(y.*sin(ALPHA.*y).*2.0+ALPHA.*s.^2.*cos(ALPHA.*y)).*2.0+A.*ALPHA.*s.*exp(ALPHA.^2.*s.^2.*(-1.0./2.0)).*(y.*sin(ALPHA.*y).*2.0+ALPHA.*s.^2.*cos(ALPHA.*y)).*4.0-A.^2.*ALPHA.^2.*s.*exp(-ALPHA.^2.*s.^2).*(exp(-ALPHA.^2.*s.^2).*cos(ALPHA.*y.*2.0)-1.0)+A.*ALPHA.^2.*s.^3.*exp(ALPHA.^2.*s.^2.*(-1.0./2.0)).*cos(ALPHA.*y).*4.0-A.^2.*ALPHA.^2.*s.*exp(-ALPHA.^2.*s.^2).*cos(ALPHA.*y.*2.0).*(exp(-ALPHA.^2.*s.^2)-1.0);

    D2Q = sum(fi(sigma, Y, A, ALPHA));
    dD2Q_ds = sum(dfi_ds(sigma, Y, A, ALPHA));
end