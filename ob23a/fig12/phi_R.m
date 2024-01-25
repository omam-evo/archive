function[phi_R, D_Q] = phi_R(MU, LAM, N, A, ALPHA, sigma, R)
    sigma_Y = R/sqrt(N);
    
    D2_Q = 2*N*sigma.^4 + 4*R.^2.*sigma.^2 + ...
        N*A^2/2 * (1-exp(-(ALPHA*sigma).^2)) .* (1-exp(-ALPHA^2*(sigma.^2 + 2*sigma_Y.^2))) + ...
        2*N*A*ALPHA^2*sigma.^2 .* exp(-0.5*ALPHA^2*(sigma.^2+sigma_Y.^2)) .* (sigma.^2 + 2*sigma_Y.^2);

    D_Q = sqrt(D2_Q);
    phi_R = e_vartheta_a_b(MU/LAM, 1, 0) * sigma.^2 .* (2*R.^2) ./ D_Q  .* (2 + ALPHA^2*A*exp(-0.5*ALPHA^2*(sigma.^2+R.^2/N))) - N*sigma.^2/MU;

end