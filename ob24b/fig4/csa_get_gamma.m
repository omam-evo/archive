function [gamma] = csa_get_gamma(THETA,C,D,N)
    cv = e_vartheta_a_b(THETA,1,0);    
    a = sqrt(2*N)*cv;
    b = (1-C)/D/(C+(1-C)*a/N);
    % b_tilde = 1/(1+sqrt(2)*cv); %large N, C=1/sqrt(N), D=sqrt(N)
    gamma = sqrt(0.5*(sqrt(1+b^2)-b+1));
end

