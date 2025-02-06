function [y, f_g, r_g, sigma_g, gen, y_g] = muComLam_sSA(SET_META_EP, FIT, N, MU, LAM, Y_0, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, G, verbose)

    %% ES
    f_g = nan*zeros(G,1);
    sigma_g = nan*ones(G,1);
    r_g = nan*ones(G,1);
    y_g = nan*ones(N, G);
  
    y = Y_0;
    s_rec = SIGMA_0;
    f_g(1) = FIT.get_f(y);
    r_g(1) = norm(y);
    sigma_g(1) = s_rec;
    y_g(:,1) = y;
    
    for g = 2:G
        
        %% Vectorized
        if SET_META_EP==0
            s_lam = repmat(s_rec * exp(TAU * randn(1, LAM)), N, 1); % sigma indep. of dim
        else
            s_lam = repmat(s_rec * (1 + TAU * randn(1, LAM)), N, 1);
        end
        z_n_lam = randn(N, LAM);
        y_n_lam = repmat(y, 1, LAM) + s_lam.*z_n_lam;
        F_lam = FIT.get_f(y_n_lam, g, norm(y-FIT.y_hat));
 
        %% Selection
        [~, idx] = sort(F_lam, 'ascend');
        z_n_mu = z_n_lam(:, idx(1:MU));
        y_n_mu = y_n_lam(:, idx(1:MU));
        s_mu = s_lam(1, idx(1:MU)); 
        
        %% Recombine & Save
        y = mean(y_n_mu, 2);
        s_rec = mean(s_mu, 2); 
        
        r_g(g) = norm(y - FIT.y_hat);
        sigma_g(g) = s_rec;
        f_g(g)= FIT.get_f(y);
        y_g(:,g) = y;
  
        %% Check
        if s_rec < SIGMA_STOP
            fprintf('\t\t SIGMA_STOP. Gen: %i, f: %d, r: %d, sig: %d \n', g-1, f_g(g), r_g(g), sigma_g(g))
            break
        end
        if r_g(g) < R_STOP || f_g(g) < F_STOP
            fprintf('\t\t F_R_STOP. Gen: %i, f: %d, r: %d, sig: %d \n', g-1, f_g(g), r_g(g), sigma_g(g))
            break
        end
        if verbose == 1
            fprintf('Gen: %i, f: %d, r: %d, sig: %d \n', g-1, f_g(g), r_g(g), sigma_g(g))
        end

    end
    if g == G
        fprintf('\t\t N_G limit reached. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g, f_g(g), r_g(g), sigma_g(g))
    end 
    gen = g;
end