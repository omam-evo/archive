function [y_n_rec, F_opt_vec, delta_r_vec, sigma_vec, gen] = muComLam_sSA(fit, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, R_STOP, N_G, verbose)
    
    % Values over generations
    F_opt_vec = nan*zeros(N_G+1,1);
    sigma_vec = nan*ones(N_G+1,1);
    delta_r_vec = nan*ones(N_G+1,1);

    % Initialize
    if length(Y_0) == 1
        y_n_mu = Y_0 * ones(N, 1);
    else
        y_n_mu = Y_0;
    end
    y_n_mu = repmat(y_n_mu, 1, MU);
    s_mu = SIGMA_0 * ones(1, MU);

    for g = 1:N_G+1

        % Recombine
        y_n_rec = mean(y_n_mu, 2);
        s_rec = mean(s_mu, 2);

        % Evaluate
        delta_r_vec(g) = norm(y_n_rec - Y_HAT);
        sigma_vec(g) = s_rec;
        F_opt_vec(g)= fit(y_n_rec);

        if s_rec < SIGMA_STOP
            fprintf('\t\t SIGMA_STOP. Gen: %i, f: %d, r: %d, sig: %d \n', g-1, F_opt_vec(g), delta_r_vec(g), sigma_vec(g))
            break
        end
        if delta_r_vec(g) < R_STOP
            fprintf('\t\t R_STOP. Gen: %i, f: %d, r: %d, sig: %d \n', g-1, F_opt_vec(g), delta_r_vec(g), sigma_vec(g))
            break
        end
        if verbose == 1
            fprintf('Gen: %i, f: %d, r: %d, sig: %d \n', g-1, F_opt_vec(g), delta_r_vec(g), sigma_vec(g))
        end
      
      %% Sequential
%       y_n_lam = nan*ones(N_DIM, LAM); s_lam = nan*ones(1, LAM); F_lam = nan*ones(1, LAM);
%       for L = 1:LAM
%         s_lam(1,L) = s_rec * exp(TAU * randn(1));
%         y_n_lam(:, L) = y_n_rec + s_lam(1,L) * randn(N_DIM,1);
%         F_lam(1,L) = fit(y_n_lam(:, L));
%       end

        %%   Vectorized
        s_lam = repmat(s_rec * exp(TAU * randn(1, LAM)), N, 1); % sigma indep. of dim
        y_n_lam = repmat(y_n_rec, 1, LAM) + s_lam.*randn(N, LAM);
        F_lam = fit(y_n_lam);

        %% Selection
        [~, idx] = sort(F_lam, 'ascend');
        y_n_mu = y_n_lam(:, idx(1:MU));
        s_mu = s_lam(1, idx(1:MU));


    end
    if g == N_G
        fprintf('N_G limit reached. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g-1, F_opt_vec(g), delta_r_vec(g), sigma_vec(g))
    end 
    gen = g-1;
end