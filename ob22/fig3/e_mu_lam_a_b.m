
%% TES Eq. 5.112, generalized progress coefficient for (mu/mu,lam)
% > e_mu_lam_a_b(500, 1000, 1, 0) ==> working
% > e_mu_lam_a_b(600, 1200, 1, 0) ==> NOT working; nchoosek(1200, 600) = INF

function e = e_mu_lam_a_b(MU, LAM, a, b)
   
    %% Prefactor (pop. dependent)
    num_sampling_points = 10000;
    c = (LAM-MU)/(2*pi)^((a+1)/2)*nchoosek(LAM, MU);
    
    %% Integrand 1
    % Population independent; normal distribution & std.normal for a=1 (c_mu_lam)
    % x_lim_I1: boundary of relevant integration range independent of population; integrand zero outside
    % fac1 = 6: 6 standard deviations as range; term x^b also suppressed by exp(-x^2)
    I1 = @(x) x.^b .* exp(-x.^2*(a+1)/2);
    x_0 = 0; 
    sigma = 1/sqrt((a+1)/2);
    fac1 = 6;
    x_lim_I1 = [-fac1*sigma + x_0, fac1*sigma + x_0];
    x = linspace(x_lim_I1(1), x_lim_I1(end), num_sampling_points);
    
    %% Integrand 2
    % Population dependent; 
    %  > may be sharp for mu/lam=0.5 and large pop.
    %  > or heavily skewed for mu/lam ~ 0 or ~1
    I2 = @(x) c * normcdf(x).^(LAM-MU-1) .* (1-normcdf(x)).^(MU-a);
    I = I2(x);
      
    % Find maximum value (locate center of peak, max be at boundary of x, i.e. at -4*sigma, or 4*sigma)
    % If at boundary: no worries bc. I1 is zero anyway
    [I_max, id_max] = max(I2(x)); 
    % Find closest id where integrand value dropped to 1/2 of maximum => peak width measure
    [~, id_sort] = sort(abs(I_max/2-I));
    % Find delta_x (half width of half maximum)
    delta_x = abs(x(id_sort(1))-x(id_max));
    % Define relevant I2 range using delta_x and some prefactor fac2 (e.g. 5)
    % I2 is zero outside, or not zero if heavily skewed, e.g. (1,300), but other itegrand will be zero
    fac2 = 5;
    x_lim_I2 = [x(id_max) - fac2*delta_x, x(id_max) + fac2*delta_x];
    
    % Bugfix needed for (1,2)-ES, I2 is const function and range of I1 is relevant
    if x_lim_I2(1) == x_lim_I2(2)
        x_lim_I2 = x_lim_I1;
    end
    
    %% Define integration range of c*I1*I2
    % max(smallest values), min(largest values) => picking only relevant range such that both integrands approx. not zero
    x_lim = [max([x_lim_I1(1),x_lim_I2(1)]), min([x_lim_I1(2),x_lim_I2(2)])]; 
    
    %% Plot integrands
%     figure();
%     x = linspace(x_lim(1), x_lim(end), num_sampling_points);
%         % I1
%         subplot(3,1,1); hold on;  ylabel('I1')
%         plot(x, I1(x), 'r');
%         % I2
%         subplot(3,1,2); hold on;  ylabel('I2')
%         plot(x, I2(x), 'r');
%         % Product
%         subplot(3,1,3); hold on;  ylabel('I1*I2')
%         plot(x, I1(x).*I2(x));
    
    %% Numeric Integration
    syms x
    e = double(vpaintegral(x^b * exp(-x^2*(a+1)/2) * c * normcdf(x)^(LAM-MU-1) * (1-normcdf(x))^(MU-a), x, x_lim));
    
end