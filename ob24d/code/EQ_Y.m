function [EQ,dEQ_ds] = EQ_Y(A, ALPHA, sigma, Y)

    fi = @(s,y,A,ALPHA) s.^2 + A * cos(ALPHA*y).*(1 - exp(-0.5*ALPHA^2*s.^2));
    dfi_ds = @(s,y,A,ALPHA) 2*s + A*ALPHA^2*s.*exp(-0.5*ALPHA^2*s.^2).*cos(ALPHA*y);

    EQ = sum(fi(sigma, Y, A, ALPHA));
    dEQ_ds = sum(dfi_ds(sigma, Y, A, ALPHA));
end