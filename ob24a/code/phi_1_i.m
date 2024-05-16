function[phi1] = phi_1_i(MU, LAM, A, ALPHA, sigma, Y, OPTIONS)
    if any(strcmpi(OPTIONS,'hybrid'))
        [D2Q, ~] = D2Q_R(A, ALPHA, sigma, norm(Y), length(Y));
    else
        [D2Q, ~] = D2Q_Y(A, ALPHA, sigma, Y);     
    end
    if any(strcmpi(OPTIONS,'set_exp_zero'))
        bool_exp = 0;
    else
        bool_exp = 1;
    end 
    D_Q = sqrt(D2Q);
    exp_di = bool_exp*exp(-(ALPHA*sigma).^2/2)*ALPHA*A.*sin(ALPHA*Y);
    phi1 = e_vartheta_a_b(MU/LAM, 1, 0) * sigma.^2 ./ D_Q .* (2*Y + exp_di);
end
        