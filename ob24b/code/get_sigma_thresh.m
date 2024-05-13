function [sigma_esc] = get_sigma_thresh(ALPHA, A)
    %% SOL BY HGB
    % x = tan(x) => first non-trivial solution
    x0 = 4.49340945790906;
    sigma_esc = 1/ALPHA*sqrt(2*log(-ALPHA^2*A*sin(x0)/(2*x0)));

    %% FIRST SEMI-ANALYTIC APPROACH
    %     G = @(y) y.*sin(ALPHA*y);
    %     dG = @(y) sin(ALPHA*y) + ALPHA*y.*cos(ALPHA*y);
    %     init = 3/4;
    %     y0 = fzero(dG, init);
    %     f0 = G(y0);
    %     if f0>=0
    %         error('Evaluated y0 yields pos gain (should be neg.)');
    %     else
    %         f0 = G(y0);
    %     end
    %     sigma_esc = sqrt(2)/ALPHA*sqrt(log((-f0)*ALPHA/2/y0^2*A));
    
end