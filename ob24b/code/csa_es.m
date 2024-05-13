%% CSA to call externally
% C = 1/sqrt(N_DIM);
% D = sqrt(N_DIM);
% E_chi = sqrt(N_DIM)*(1 - 1/(4*N_DIM) + 1/(21*N_DIM^2));

function [y, f_g, r_g, sigma_g, gen, y_g] = csa_es(FIT, N, MU, LAM, Y_0, SIGMA_0, S_PATH_0, C, D, E_chi, SIGMA_STOP, R_STOP, F_STOP, G, verbose)
    
    f_g = nan*zeros(G,1);
    sigma_g = nan*ones(G,1);
    r_g = nan*ones(G,1);
    y_g = nan*ones(N,G); 

    y = Y_0;
    s = S_PATH_0 * ones(N,1);
    sigma = SIGMA_0;
    
    % --- needed for non-vectorized version ---
    %z_tilde = nan*ones(N, LAM); 
    %y_tilde = nan*ones(N, LAM); 
    %f_tilde = nan*ones(1, LAM);
    
    for g = 1:G
        
        z_tilde = randn(N, LAM);
        y_tilde = repmat(y, 1, LAM) + sigma*z_tilde;
        f_tilde = FIT.get_f(y_tilde, g, norm(y-FIT.y_hat));
        
%         for i=1:LAM
%             z_tilde(:, i) = randn(N, 1);
%             y_tilde(:, i) = y_p + sigma*z_tilde(:, i);
%             f_tilde(i) = FIT(y_tilde(:,i));
%         end

        [~, idx] =  sort(f_tilde, 'ascend');
        y_mu = y + sigma*z_tilde(:, idx(1:MU));
        y = mean(y_mu, 2);
        z_p = mean(z_tilde(:, idx(1:MU)), 2);
        %y_p = y_p + sigma*z_p; % old
        
        %% Update sigma
        s = (1-C)*s + sqrt(MU*C*(2-C))*z_p;
        %sigma = sigma * exp((norm(s)-E_chi) / (D*E_chi)); % Hansen
        sigma = sigma * exp((norm(s)^2-E_chi) / (2*D*E_chi)); %HGB, AR
        
        f_g(g) = FIT.get_f(y);
        sigma_g(g) = sigma;
        r_g(g) = norm(y-FIT.y_hat);
        y_g(:,g) = y;
        
        if verbose == 1
            fprintf('Gen: %i, fit: %d, sigma: %d \n', g-1, f_g(g), sigma_g(g));
        end
        if sigma_g(g) < SIGMA_STOP
            fprintf('\t\t SIGMA_STOP. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g-1, f_g(g), r_g(g), sigma_g(g))
            break
        end
        if r_g(g) < R_STOP || f_g(g) < F_STOP
            fprintf('\t\t F_R_STOP. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g-1, f_g(g), r_g(g), sigma_g(g))
            break
        end
         
    end
    if g == G
        fprintf('\t\t N_G limit reached. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g, f_g(g), r_g(g), sigma_g(g))
    end 
    gen = g;
end
    