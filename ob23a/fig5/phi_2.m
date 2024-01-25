function[phi_2] = phi_2(MU, LAM, A, ALPHA, sigma, Y, e_11, e_20)
    bool_exp = 1;
    [phi, D_Q] = phi_1(MU, LAM, A, ALPHA, sigma, Y, bool_exp);

    %% Get Phi2
    phi_2 = 2*Y.*phi ...
        - sigma^2 / MU * (1 +(2*Y).^2 * sigma^2 ./ D_Q^2 * (e_11 + (MU-1)*e_20));
end