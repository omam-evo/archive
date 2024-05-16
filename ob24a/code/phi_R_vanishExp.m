function[phi_R, D_Q] = phi_R_vanishExp(MU, LAM, N, A, ALPHA, sigma, R)
    ZERO = 0;
    phi_R = e_vartheta_a_b(MU/LAM, 1, 0) * sigma.^2 .* (2*R.^2) ./ sqrt(D2Q_R_terms(A, ALPHA, sigma, R, N))  .* ...
        (2 + ALPHA^2*A*exp(-0.5*ALPHA^2*(sigma.^2+R.^2/N))) - N*sigma.^2/MU;
end