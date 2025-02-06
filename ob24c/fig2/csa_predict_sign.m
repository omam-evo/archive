function [res] = csa_predict_sign(MU,LAM,N,C,D,SET_B_SIMPLI)
    %% predict sign*
    cv = e_vartheta_a_b(MU/LAM,1,0);
    a = sqrt(2*N)*cv;
    if SET_B_SIMPLI == 0
        b = (1-C)/D/(C+(1-C)*a/N);
    else
        b = (1-C)/(C*D);
    end
    
    %% numeric
    % phinorm = @(x) a-x^2/(2*MU);
    % s = @(x) a*MU/x^2*(a-x^2/MU);
    % f = @(x) phinorm(x) + b*s(x);
    % sign_num = fzero(f, SIGNZERO);
    
    % sign_analytic_ab = sqrt(a*MU*(sqrt(1+b^2)-b+1));
    res = (2*N)^(1/4)*(cv*MU)^(1/2)*sqrt(sqrt(1+b^2)-b+1);
end

