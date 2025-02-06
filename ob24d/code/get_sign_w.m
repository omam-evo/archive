function [sign0_w] = get_sign_w(N, A, ALPHA, MU, LAM)
    sign0_w = (512*N)^(1/4)*sqrt(e_vartheta_a_b(MU/LAM,1,0)*MU)/ALPHA/sqrt(A)*sqrt(lambertw(0, ALPHA^3*A^(3/2)/8));
end