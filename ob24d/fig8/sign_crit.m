function [sign0_phi, sign0_w, sign0_polyn1] = sign_crit(MU,LAM,N,A,ALPHA)
%     A=2;
%     ALPHA=2*pi;
%     N=100;
%     MU=100;
%     LAM=200;
    
    e_10 = e_vartheta_a_b(MU/LAM, 1, 0);
    sigma_eps = sqrt(N/2)*A;
    R = sqrt(sigma_eps*N/(4*e_10*MU));
    ZERO = 0;

    x0 = 0.5*sqrt( sqrt(N*(8*e_10^2*MU^2 + N)) - N);
    f1 = @(x) 1+0.5*ALPHA^2*A*exp(-0.5*(ALPHA*x)^2*A/(sqrt(32*N)*e_10*MU))-sqrt(1+x^2/(4*e_10^2*MU^2)+x^4/(8*N*e_10^2*MU^2));
    f2 = @(x) 0.5*ALPHA^2*A*exp(-0.5*(ALPHA*x)^2*A/(sqrt(32*N)*e_10*MU))-x^4/(16*N*e_10^2*MU^2);
    sign0_f1 = fzero(f1,20);
    sign0_f2 = fzero(f2,20);

    D2_Q = @(s) 4*R.^2.*s.^2 + 2*N*s.^4 + ...
            N*A^2/2 * (1-exp(-(ALPHA*s).^2)) .* (1-exp(-ALPHA^2*(s.^2 + 2*(R/sqrt(N)).^2))) + ...
            2*N*A*ALPHA^2*s.^2 .* exp(-0.5*ALPHA^2*(s.^2+(R/sqrt(N)).^2)) .* (s.^2 + 2*(R/sqrt(N)).^2);

    phi_R = @(s) e_10 * (s.*R/N).^2 .* (2*R.^2) ./ sqrt(D2_Q(s.*R/N))  .* ...
            (2 + ALPHA^2*A*exp(-0.5*(ALPHA*R)^2*((s/N).^2+1/N))) - N*(s.*R/N).^2/MU;

    sign0_phi = fzero(phi_R,x0);

    %% solution a*exp(-b*x^2) = x^4
    syms a b x
    eqn = a*exp(-b*x^2)==x^4;
    res = solve(eqn,x, 'ReturnConditions',true);
    a = 8*ALPHA^2*A*N*e_10^2*MU^2;
    b = 0.5*ALPHA^2*A/(sqrt(32*N)*e_10*MU);

    k = 0;
    s1 = double(subs(res.x(1)));
    s2 = double(subs(res.x(2)));
    s3 = double(subs(res.x(3)));
    s4 = double(subs(res.x(4)));
    
    
    %sign0_lam = (2^(1/2)*lambertw(k, (a^(1/2)*b)/2)^(1/2))/b^(1/2);
    %sign0_w = (512*N)^(1/4)*sqrt(e_10*MU)/ALPHA/sqrt(A)*sqrt(lambertw(0, ALPHA^3*A^(3/2)/8));
    sign0_w = get_sign_w(N, A, ALPHA, MU, LAM);
    
     %% solution a*exp(-b*x^2) = x^2
%     sign0_w = (32*N)^(1/4)*sqrt(2*cv*MU)/ALPHA/sqrt(A)*sqrt(lambertw(0, 2*cv*MU*ALPHA^4*A^2/(32*N)));

    %% solution of polynomial instead of exp(..)
%     eqn = a*(1 - 1*b*x^2 + ZERO*b^2*x^4/2 -ZERO*b^3*x^6/6 + ZERO*b^4*x^8/24) == x^4 + ZERO*2*N*x^2;
%     res2 = vpasolve(eqn,x,x0);

    sign0_polyn1 = sqrt(sqrt(a)*sqrt(a*b^2 + 4) - a*b)/sqrt(2);  %x^4 + ab*x^2 - a = 0;
