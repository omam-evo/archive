%% Rework after comparing with HGB coefficient solver

function e = e_mu_lam_a_b_v2(mu, lam, a, b)
   
    %% Prefactor (pop. dependent)
    num_sampling_points = 10000;
    
    c = 1/(2*pi)^((a+1)/2);
    logfac = gammaln(lam+1)-gammaln(mu+1)-gammaln(lam-mu);
    
    I = @(x) e_integrand(x, mu, lam, a, b, logfac);
    f1 = @(x) -(a+1)*x + (lam-mu-1)*normpdf(x)/normcdf(x) - (mu-a)*normpdf(x)/(1-normcdf(x));
    x0 = fzero(f1, 0);
%     f2 = @(x) 1+x*(-(a+1)*x + (lam-mu-1)*normpdf(x)/normcdf(x) - (mu-a)*normpdf(x)/(1-normcdf(x)));
%     x0 = fzero(f2, x0_f1);
    
    x_min = x0-6;
    x_max = x0+6;
    
    e = c*integral(I, x_min, x_max);
    
    %% Plot integrands
%     x = linspace(x_min, x_max, num_sampling_points);
%     num = trapz(x,c*I(x));
%     figure(); hold on; title(['e=',num2str(e),', num=',num2str(num)]);
%         plot(x, c*I(x), 'r');    

end