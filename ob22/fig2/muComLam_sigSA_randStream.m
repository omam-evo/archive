function [y_n_rec, F_opt_vec, delta_r_vec, sigma_vec, gen] = muComLam_sigSA_randStream(rstream, fit, N_DIM, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, F_STOP, N_G, verbose)
    
    % Values over generations
    F_opt_vec = nan*zeros(N_G+1,1);
    sigma_vec = nan*ones(N_G+1,1);
    delta_r_vec = nan*ones(N_G+1,1);

    % Initialize
    if length(Y_0) == 1
        y_n_mu = Y_0 * ones(N_DIM, 1);
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
            fprintf('\t\t SIGMA_STOP. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g-1, F_opt_vec(g), delta_r_vec(g), sigma_vec(g))
            break
        end
        if F_opt_vec(g) < F_STOP
            fprintf('\t\t F_STOP. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g-1, F_opt_vec(g), delta_r_vec(g), sigma_vec(g))
            break
        end
        if verbose == 1
            fprintf('Gen: %i, fit: %d, delta_r: %d, Sigma %d \n', g-1, F_opt_vec(g), delta_r_vec(g), sigma_vec(g))
        end
      
      % Initialize storage
%       y_n_lam = nan*ones(N_DIM, LAM);
%       s_lam = nan*ones(1, LAM);
%       F_p_lam = nan*ones(1, LAM);

%       for L = 1:LAM
%         s_l = s_rec * exp(TAU * randn(1));
%         y_n = y_n_rec + s_l * randn(N_DIM,1);
%         F_lam = fit(y_n);
% 
%         s_lam(1,L) = s_l;
%         y_n_lam(:, L) = y_n;
%         F_p_lam(1,L) = F_lam;
%       end

        %%   Vectorized
        s_lam = repmat(s_rec * exp(TAU * randn(rstream, 1, LAM)), N_DIM, 1); % sigma indep. of dim
        y_n_lam = repmat(y_n_rec, 1, LAM) + s_lam.*randn(rstream, N_DIM, LAM);
        F_lam = fit(y_n_lam);

        %% Selection
        [~, idx] = sort(F_lam, 'ascend');
        y_n_mu = y_n_lam(:, idx(1:MU));
        s_mu = s_lam(1, idx(1:MU));


    end
    if N_G~=1 && g == N_G+1
        fprintf('N_G limit reached. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g-1, F_opt_vec(g), delta_r_vec(g), sigma_vec(g))
    end 
    gen = g-1;
end