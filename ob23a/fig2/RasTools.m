classdef RasTools
   properties
    f; f_d1; f_d2; f_d3; q_i;
    A; alpha;
    Y_0; 
    a_quad;
    Q_diag_quad;
    fac_sigma2norm;
    fac_q2norm;
    f_cos1; f_cos2; f_sin2; f_x2cos1; f_x1sin1; f_x1cos1; f_x3;
    EV_cos1; EV_cos2; EV_sin2; EV_x2cos1; EV_x1sin1; EV_x1cos1; EV_x3; VAR_cos1; VAR_sin1;

   end
   methods
       
       function [obj] = RasTools()
            %% Derivatives
            obj.f =  @(x, alpha, a) sum(a - a * cos(alpha * x) + x.^2, 1);
            obj.q_i = @(a, alpha, y, x) x.^2 +2*x*y - a*cos(alpha*y)*cos(alpha*x) + a*sin(alpha*y)*sin(alpha*x) + a*cos(alpha*y);
            obj.f_d1 = @(x, alpha, a) 2*x + a * alpha * sin(alpha*x);
            obj.f_d2 = @(x, alpha, a) 2 + a * (alpha)^2 * cos(alpha*x);
            obj.f_d3 = @(x, alpha, a) -a * (alpha)^3 * sin(alpha*x);
            
            %% Functions and expected values
            obj.f_cos1 = @(a,x) cos(a*x);
            obj.EV_cos1 = @(a,mu,s) exp(-(a*s)^2/2)*cos(a*mu);
            
            obj.f_cos2 = @(a,x) cos(a*x).^2;
            obj.EV_cos2 = @(a,mu,s) 1/2*(1 + obj.EV_cos1(2*a,mu,s));
            
            obj.f_sin2 = @(a,x) sin(a*x).^2;
            obj.EV_sin2 = @(a,mu,s) 1/2*(1 - obj.EV_cos1(2*a,mu,s));
            
            obj.f_x2cos1 = @(a,x) x.^2.*cos(a*x);
            obj.EV_x2cos1 = @(a,mu,s) exp(-(a*s)^2/2)*(mu^2*cos(a*mu) + s^2*cos(a*mu) - 2*a*mu*s^2*sin(a*mu) - a^2*s^4*cos(a*mu));
            
            obj.f_x1sin1 = @(a,x) x.*sin(a*x);
            obj.EV_x1sin1 = @(a,mu,s) exp(-(a*s)^2/2)*(a*s^2*cos(a*mu)+mu*sin(a*mu));

            obj.f_x1cos1 = @(a,x) x.*cos(a*x);
            obj.EV_x1cos1 = @(a,mu,s) exp(-(a*s)^2/2)*(mu*cos(a*mu) - a*s^2*sin(a*mu));
            
            obj.VAR_cos1 = @(a,mu,s) obj.EV_cos2(a,mu,s) - obj.EV_cos1(a,mu,s)^2;
            obj.VAR_sin1 = @(a,mu,s) obj.EV_sin2(a,mu,s) - 0;
            
            obj.f_x3 = @(x) x.^3;
            obj.EV_x3 = @(mu,s) mu^3 + 3*mu*s^2;
            
       end
       function[T] = get_taylor(obj, alpha, A, order, x_0)
           syms x
           T(x) = taylor(obj.f(x, alpha, A), x, 'ExpansionPoint', x_0, 'Order', order+1);
       end
      function [x_max, x_min] = get_extrema(obj, alpha, A, order_n)
        options = optimset('Display', 'off');
        eps = 1e-3;
        f_d1_A = @(x) obj.f_d1(x,alpha,A);
          
        f_x_extr_est = @(x_grid,alpha,A) (x_grid*(alpha)^2*A*cos(alpha*x_grid) - alpha*A*sin(alpha*x_grid)) / (2 + (alpha)^2*A*cos(alpha*x_grid));
          
        x_max = nan*ones(length(order_n), 1);
        x_min = nan*ones(length(order_n), 1);
        
        x_threshold = pi*A;
        order_n_unique = unique(order_n);
        
        for i = 1:length(order_n_unique)
            ord = order_n_unique(i);
            
            %% MAX
            x_grid = 1*ord-0.5;
            if ord ~= 0 && x_grid <= ceil(x_threshold)
                x_extr_est = f_x_extr_est(x_grid, alpha, A);
                [xs, fs] = fsolve(f_d1_A, x_extr_est, options);
                if abs(fs) < eps
                    x_max(order_n==ord) = xs;
                end
            else
                x_max(order_n==ord) = nan;
            end
            
            %% MIN
            x_grid = 1*ord;
            if x_grid <= ceil(x_threshold)
                x_extr_est = f_x_extr_est(x_grid,alpha,A);
                [xs, fs] = fsolve(f_d1_A, x_extr_est, options);
                if abs(fs) < eps
                    x_min(order_n==ord) = xs;
                end
             else
                x_min(order_n==ord) = nan;
            end   
        end
        
      end
   
      function [order_y_min, order_y_max] = get_all_extrema(obj, alpha, A, N_DIM)
            order_list = [0,1,2,3,4,5,6,7];
            order_y_min = nan*ones(1,length(order_list));
            order_y_max = nan*ones(1,length(order_list));

            for i = 1:length(order_list)
                order_dim = zeros(N_DIM, 1); 
                order_dim(1) = order_list(i);
                [Y_0_max, Y_0_min] = obj.get_extrema(alpha, A, order_dim);
                order_y_min(i) = Y_0_min(1);
                order_y_max(i) = Y_0_max(1);
            end
            order_y_min = order_y_min(~isnan(order_y_min));
            order_y_max = order_y_max(~isnan(order_y_max));
      end
        
    %% Exp-Value of Rastrigin Q-Gain
    function [E_q] = expval_q_gain(obj, A, alpha, Y, sigma)
        exp_cos_x = obj.EV_cos1(alpha, 0, sigma);
        cos_y = cos(alpha*Y);
        E_q = sigma^2 + A * cos_y - A * cos_y * exp_cos_x;
    end
    
    %% E[Q|x_i] linear
    function [E_q_fix] = expval_q_gain_fix(obj, A, alpha, x_i, y_i, Y_j_not_i, sigma, order)
        T = x_i*obj.f_d1(y_i, alpha, A);
        if order == 2
            T = T + x_i^2/2*obj.f_d2(y_i, alpha, A);
        end
        % T = T + x_i^3/6*obj.f_d3(y_i, alpha, A);
        % T = x_i*2*y_i;
        E_q_fix = T + sum(expval_q_gain(obj, A, alpha, Y_j_not_i, sigma));
    end
    %% E[Q|x_i] non-linear
    function [E_q_fix] = expval_q_gain_fix_nonlin(obj, A, alpha, x_i, y_i, Y_j_not_i, sigma)
        E_q_fix = obj.q_i(A, alpha, y_i, x_i) + sum(expval_q_gain(obj, A, alpha, Y_j_not_i, sigma));
    end
    %% E[Q|x_i] test
    function [E_q_fix] = expval_q_gain_fix_testing(obj, A, alpha, x_i, y_i, Y_j_not_i, sigma)
        E_q_fix = 0*x_i.^2 + 2*x_i*y_i + 0*A*cos(alpha*y_i) - 0*A*cos(alpha*y_i)*cos(alpha*x_i) + A*sin(alpha*y_i)*sin(alpha*x_i) + ...
            sum(expval_q_gain(obj, A, alpha, Y_j_not_i, sigma));
    end
    
    %% Variance: first step using Var[.] = E[(.)^2]-E[.]^2
    function [Var_q] = variance_q_gain_v1(obj, A, alpha, y_mu, y_sigma, sigma)
        
        % Initialize
        mutation_mu = 0; % mutation is centered around 0
        N_DIM = length(y_mu);

        EV_x_c1 = obj.EV_cos1(alpha, mutation_mu, sigma);
        EV_x_c2 = obj.EV_cos2(alpha, mutation_mu, sigma);
        EV_x_s2 = obj.EV_sin2(alpha, mutation_mu, sigma);
        EV_x_x1s1 = obj.EV_x1sin1(alpha, mutation_mu, sigma);    
        EV_x_x2c1 = obj.EV_x2cos1(alpha, mutation_mu, sigma);

        % Given sigma: N_DIM components
        q_i_exp = nan*ones(N_DIM, 1);
        q_i_sq_exp = nan*ones(N_DIM, 1);

        for j = 1:N_DIM
            EV_y_c1 = obj.EV_cos1(alpha, y_mu(j), y_sigma);
            EV_y_c2 = obj.EV_cos2(alpha, y_mu(j), y_sigma);
            EV_y_s2 = obj.EV_sin2(alpha, y_mu(j), y_sigma);
            EV_y_x1s1 = obj.EV_x1sin1(alpha, y_mu(j), y_sigma);

            q_i_exp(j) = expval_q_gain(obj, A, alpha, y_mu(j), sigma);
            
            q_i_sq_exp(j) = ...
                3 * sigma^4 + ... 
                4 * sigma^2 * (y_mu(j)^2+y_sigma^2) + ... 
                2 * sigma^2 * (A * EV_y_c1) + ...
                A^2 * EV_y_c2 + ...
                A^2 * EV_y_c2 * EV_x_c2 + ...
                A^2 * EV_y_s2 * EV_x_s2 + ...
                (-2 * A * EV_y_c1 * EV_x_x2c1) + ...
                (-2 * A^2 * EV_y_c2 * EV_x_c1) + ...
                4 * A * EV_y_x1s1 * EV_x_x1s1;
        end
        
        Var_q = q_i_sq_exp - q_i_exp.^2;
    end
    
    %% Variance: 2nd simplification step
    function [Var_q] = variance_q_gain_v2(obj, A, alpha, Y, sigma)
        
        ev_cos1 = obj.EV_cos1(alpha, 0, sigma);
        ev_x1sin1 = obj.EV_x1sin1(alpha, 0, sigma);
        ev_x2cos1 = obj.EV_x2cos1(alpha, 0, sigma);
        ev_sin2 = obj.EV_sin2(alpha, 0, sigma);
        
        var_cos = obj.EV_cos2(alpha, 0, sigma) - ev_cos1^2;
        var_sin = ev_sin2 - 0;
        
        Var_q = nan*ones(length(Y), 1);
        for i=1:length(Y)
            Var_q(i) = 2*sigma^4 + ... 
                4*Y(i)^2*sigma^2 + ...
                A^2 * sin(alpha*Y(i))^2 * var_sin + ...
                A^2 * cos(alpha*Y(i))^2 * var_cos + ...
                -2*A* cos(alpha*Y(i)) * ev_x2cos1 + ...
                2*A*sigma^2 * cos(alpha*Y(i)) * ev_cos1 + ... 
                4*A * Y(i)*sin(alpha*Y(i)) * ev_x1sin1;
        end
        
    end
    
    %% Standard simplified variance expression
    function [Var_q] = variance_q_gain_simplified(obj, A, alpha, Y, sigma)  
        expfac = exp(-(alpha*sigma)^2);
        Var_q = 2*sigma^4 + 4*Y.^2*sigma^2 + A^2/2*(1 - expfac)*(1-cos(2*alpha*Y)*expfac) + ...
            2*A*alpha*sigma^2*exp(-0.5*(alpha*sigma)^2)*(alpha*sigma^2*cos(alpha*Y) + 2*Y.*sin(alpha*Y));
    end
    
    %% Variance as function of R and sigma
    function [Var_q] = variance_q_gain_R(obj, A, alpha, R, N, sigma)
        mu_Y = 0;
        sigma_Y = R/sqrt(N);
    
        % Version 1
%         ev_cos1a = obj.EV_cos1(alpha, mu_Y, sigma_Y);
%         ev_cos2a = obj.EV_cos1(2*alpha, mu_Y, sigma_Y);
%         ev_ysinay = obj.EV_x1sin1(alpha, mu_Y, sigma_Y);
%         Var_q_v1 = N*2*sigma^4 + 4*R^2*sigma^2 + N*A^2/2*(1 - exp(-(alpha*sigma)^2)) + ...
%             N*A^2/2*exp(-(alpha*sigma)^2) * ev_cos2a *(exp(-(alpha*sigma)^2) - 1) + ...
%             N*2*A*alpha*sigma^2*exp(-(alpha*sigma)^2/2)*(alpha*sigma^2*ev_cos1a + 2*ev_ysinay);
        
        %% Version 2
        Var_q = 2*N*sigma.^4 + 4*R.^2.*sigma.^2 + ...
            N*A^2/2 * (1-exp(-(alpha*sigma).^2)) .* (1-exp(-alpha^2*(sigma.^2 + 2*sigma_Y.^2))) + ...
            2*N*A*alpha^2*sigma.^2 .* exp(-0.5*alpha^2*(sigma.^2+sigma_Y.^2)) .* (sigma.^2 + 2*sigma_Y.^2);
        
    end
    
    %% Variance as function of R and sigma*
    % > see function variance_q_gain_R for sigma
    % > variance_q_gain_R() and variance_q_gain_R_sigNorm() yield same result
    function [Var_q] = variance_q_gain_R_sigNorm(obj, A, alpha, R, N, sigma_norm)
        sigN = sigma_norm/N;
        
        Var_q = R^4*(2*N*sigN^4 + 4*sigN^2) + ... 
            R^4 * 2*N*A*alpha^2*sigN^2*(sigN^2+2/N)*exp(-0.5*(alpha*R)^2*(sigN^2+1/N)) + ...
            N*A^2/2 * (1-exp(-(alpha*R)^2*sigN^2)) * (1-exp(-(alpha*R)^2*(sigN^2+2/N)));
    end
    
    %% Yields variance D^2_k=(f_p*sigma)^2+D^2_i depending on each i; 
    % each comp. phi_i depends on variance of all N comp
    function [Var_q_Dk, f_p, k_i, d_i] = variance_q_gain_Df(obj, A, ALPHA, Y, sigma, bool_di)    
        
        N = length(Y);
        Var_q_Dk = nan*Y;
        k_i = nan*Y;
        d_i = nan*Y;
        f_p = nan*Y;
        for i=1:N
            Y_i = Y(i);
            Y_jni = Y(i ~= 1:N);
            D_i_sq = sum(obj.variance_q_gain_simplified(A, ALPHA, Y_jni, sigma));
            k_i(i) = 2*Y_i;
            if bool_di == 0
                d_i(i) = 0;
            else
                d_i(i) = ALPHA * A * sin(ALPHA*Y_i);
            end
            f_p(i) = k_i(i)+d_i(i);
            Var_q_Dk(i) = (f_p(i)*sigma)^2 + D_i_sq;
        end
        
    end

    %% Approx. D_Q with small Sigma y_i-dependent
    function [Var_q] = variance_q_gain_smallSig(obj, A, alpha, Y, sigma, toggle_2nd_order)    
        Var_q = 0*Y;
        for i=1:length(Y)
            Var_q(i) = (obj.f_d1(Y(i), alpha, A)*sigma)^2 + toggle_2nd_order*0.5*(obj.f_d2(Y(i), alpha, A)*sigma^2)^2;
        end      
    end
    
    %% Variance for large sigma
    function [Var_q] = variance_q_gain_largeSig(obj, A, alpha, Y, sigma)     
        Var_q = 0*Y;
        for i=1:length(Y)
            Var_q(i) = 2*sigma^4 + A^2/2 + 4*sigma^2*Y(i)^2;
%             warning('0*2*sigma^4 0*4*sigma^2*Y(i)');
        end     
    end
    
    %% Variance for small sigma & R-dependent
    function [D_q_sq] = variance_D2_Q_sphericalWithCorr(obj, A, ALPHA, N, sigma, R)
        D_q_sq = sigma^2 * R^2 * (2+ALPHA^2*A)^2 + ...
            sigma^4/2 * N * (2+ALPHA^2*A)^2;
    end
    
    %% Get Rastrigin noise strength
    function [sigma_eps_inf, sigma_eps_R] = get_sigma_eps(obj, A, alpha, N_DIM, R)
        sigma_eps_inf = sqrt(N_DIM*A^2/2);
        mu = 0;
        s = R/sqrt(N_DIM);
        sigma_eps_R = nan*R;
        % sigma_eps_R = sqrt(N_DIM*A^2*(1/2 + 1/2*exp(-2*(alpha*s).^2) - exp(-(alpha*s).^2)));
        for i=1:length(R)
           sigma_eps_R(i) = sqrt(N_DIM*A^2 * obj.VAR_cos1(alpha, mu, s(i))); 
        end
    end
    
        %% Expected quality gain for (1,lam) => deprec?
    function [M_q, S_q, pdf_x, pdf_y] = q_gain_v2(obj, A, alpha, sigma, Y, x_i, fixed_comp_id)
        
        N_DIM = length(Y);
        
        %% Calculate
        if ~isnan(fixed_comp_id)
            mask = ones(N_DIM,1);
            mask(fixed_comp_id) = 0;
            Y_j_not_i = Y(logical(mask));
            y_i = Y(fixed_comp_id);
            E_q = expval_q_gain_fix(obj, A, alpha, x_i, y_i, Y_j_not_i, sigma, 1);
            %% SET VAR
            Var_q = variance_q_gain_v1(obj, A, alpha, Y_j_not_i, 0, sigma);
%             Var_q = variance_q_gain_v2(obj, A, alpha, Y_j_not_i, sigma);
        else
            E_q = expval_q_gain(obj, A, alpha, Y, sigma);
            %% SET VAR
            Var_q = variance_q_gain_v1(obj, A, alpha, Y, 0, sigma);
%             Var_q = variance_q_gain_v2(obj, A, alpha, Y, sigma);
        end
        
        M_q = sum(E_q);
        S_q = sqrt(sum(Var_q));
        
        %% Prepare Plot
        pd = makedist('Normal','mu', M_q, 'sigma', S_q);
        pdf_x = linspace(pd.mu-3*pd.sigma, pd.mu+3*pd.sigma, 100);
        pdf_y = pdf(pd, pdf_x);

    end
    
    end % methods
end