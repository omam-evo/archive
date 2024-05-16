function [psi,dDQ,DQ,dEQ] = psi_Y(A, ALPHA, TAU, sigma, Y, e_10, e_11) 
            
    [~, dEQ] = EQ_Y(A, ALPHA, sigma, Y);
    [D2Q, dD2Q] = D2Q_Y(A, ALPHA, sigma, Y);

    DQ = sqrt(D2Q);
    dDQ = dD2Q/DQ/2;

    psi = TAU^2*(0.5 + e_11*sigma.*dDQ./DQ - e_10*sigma.*dEQ./DQ);  
     
end

