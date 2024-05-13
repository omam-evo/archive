%% Analogous to SMN C.23; SAR on general function with R-dependent EQ and D2Q
function [psi,dDQ,DQ,dEQ] = psi_R_metaep(A, ALPHA, TAU, sigma, R, N, e_10, e_11)

    [~, dEQ] = EQ_R(A, ALPHA, sigma, R, N);
    [D2Q, dD2Q] = D2Q_R(A, ALPHA, sigma, R, N);

    DQ = sqrt(D2Q);
    dDQ = dD2Q/DQ/2;

    psi = TAU^2*(e_11*sigma*dDQ/DQ - e_10*sigma*dEQ/DQ);  

end