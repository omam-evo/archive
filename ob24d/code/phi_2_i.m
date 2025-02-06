function[phi_2] = phi_2_i(MU, LAM, A, ALPHA, sigma, Y, e_11, e_20, OPTIONS)
    if any(strcmpi(OPTIONS,'hybrid'))
        [D2Q, ~] = D2Q_R(A, ALPHA, sigma, norm(Y), length(Y));
    else
        [D2Q, ~] = D2Q_Y(A, ALPHA, sigma, Y);     
    end
    p1 = phi_1_i(MU, LAM, A, ALPHA, sigma, Y, OPTIONS);
    phi_2 = 2*Y.*p1 - sigma.^2 / MU .* (1 +(2*Y).^2 * sigma.^2 ./ D2Q * (e_11 + (MU-1)*e_20));
end