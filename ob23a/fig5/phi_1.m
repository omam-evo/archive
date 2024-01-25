function[phi_1, D_Q] = phi_1(MU, LAM, A, ALPHA, sigma, Y, bool_exp)

    D2_Q = 2*sigma^4 + 4*Y.^2*sigma^2 + A^2/2*(1 - exp(-(ALPHA*sigma)^2))*(1-cos(2*ALPHA*Y)*exp(-(ALPHA*sigma)^2)) + ...
            2*A*ALPHA*sigma^2*exp(-0.5*(ALPHA*sigma)^2)*(ALPHA*sigma^2*cos(ALPHA*Y) + 2*Y.*sin(ALPHA*Y));
    D_Q = sqrt(sum(D2_Q));
    
    ki = 2*Y;
    exp_di = bool_exp*exp(-(ALPHA*sigma)^2/2)*ALPHA*A*sin(ALPHA*Y);
    phi_1 = e_vartheta_a_b(MU/LAM, 1, 0) * sigma^2 / D_Q * (ki + exp_di);
end