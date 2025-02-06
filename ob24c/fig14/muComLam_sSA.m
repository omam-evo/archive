function [y, f_g, r_g, sigma_g, gen, F_lam] = muComLam_sSA(fit, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, G, verbose)
    
    % Values over generations
    f_g = nan*zeros(G,1);
    sigma_g = nan*ones(G,1);
    r_g = nan*ones(G,1);

    y = Y_0;
    sigma = SIGMA_0;

    for g = 1:G

        % Evaluate
        r_g(g) = norm(y - Y_HAT);
        sigma_g(g) = sigma;
        f_g(g)= fit(y);

        if sigma < SIGMA_STOP
            fprintf('\t\t SIGMA_STOP. Gen: %i, f: %d, r: %d, sig: %d \n', g, f_g(g), r_g(g), sigma_g(g))
            break
        end
        if r_g(g) < R_STOP
            fprintf('\t\t R_STOP. Gen: %i, f: %d, r: %d, sig: %d \n', g, f_g(g), r_g(g), sigma_g(g))
            break
        end
        if verbose == 1
            fprintf('Gen: %i, f: %d, r: %d, sig: %d \n', g, f_g(g), r_g(g), sigma_g(g))
        end
      
      %% Sequential
%       y_n_lam = nan*ones(N_DIM, LAM); s_lam = nan*ones(1, LAM); F_lam = nan*ones(1, LAM);
%       for L = 1:LAM
%         s_lam(1,L) = s_rec * exp(TAU * randn(1));
%         y_n_lam(:, L) = y_n_rec + s_lam(1,L) * randn(N_DIM,1);
%         F_lam(1,L) = fit(y_n_lam(:, L));
%       end

        %%   Vectorized
        s_lam = repmat(sigma * exp(TAU * randn(1, LAM)), N, 1); % sigma indep. of dim
        y_n_lam = repmat(y, 1, LAM) + s_lam.*randn(N, LAM);
        F_lam = fit(y_n_lam);

        %% Selection
        [~, idx] = sort(F_lam, 'ascend');
        y_n_mu = y_n_lam(:, idx(1:MU));
        s_mu = s_lam(1, idx(1:MU));

        %% Recombine
        y = mean(y_n_mu, 2);
        sigma = mean(s_mu, 2);

    end
    if g == G
        fprintf('N_G limit reached. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g, f_g(g), r_g(g), sigma_g(g))
    end 
    gen = g;
end