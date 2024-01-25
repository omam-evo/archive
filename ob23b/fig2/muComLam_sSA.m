function [y, f_g, r_g, sigma_g, gen, y_g] = muComLam_sSA(FIT, A, ALPHA, N, MU, LAM, Y_0, SIGMA_0, TAU, SIGMA_STOP, R_STOP, F_STOP, G, verbose)

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
        s_lam = repmat(s_rec * exp(TAU * randn(1, LAM)), N, 1); % sigma indep. of dim
        z_n_lam = randn(N, LAM);
        y_n_lam = repmat(y, 1, LAM) + s_lam.*z_n_lam;
        F_lam = FIT.get_f(y_n_lam);
 
        %% Selection
        [~, idx] = sort(F_lam, 'ascend');
        z_n_mu = z_n_lam(:, idx(1:MU));
        y_n_mu = y_n_lam(:, idx(1:MU));
        s_mu = s_lam(1, idx(1:MU)); 
        
        %% Q GAIN HISTOGRAM
        F_y0 = FIT.get_f(Y_0);
        EQ = EQ_Y(A, ALPHA, SIGMA_0, Y_0);
        D2Q = D2Q_Y(A, ALPHA, SIGMA_0, Y_0);
        DQ = sqrt(D2Q);
        
        figure; hold on;
        xlabel('$f(\mathbf{y}+\mathbf{x}_k)-f(\mathbf{y})$');
        ylabel('PDF');
        histogram(F_lam-F_y0, 'normalization', 'pdf');
        x = linspace(EQ-3*DQ, EQ+3*DQ, 1001);
        pdf_x = pdf('Normal',x,EQ,DQ);
        plot(x, pdf_x, 'r-');
%         xline(EQ, 'r-');
%         xline(EQ+DQ, 'r--');
%         xline(EQ-DQ, 'r--');
        dummy = 1;
        myfigstyle(gcf,9,5,9,9)
        
        %% Hypothesis test
        dist = makedist('normal','mu',EQ,'sigma',DQ);
        [h,p] = adtest(F_lam-F_y0,'Distribution', dist);
        if h==0
            fprintf(' fail to reject H0 (normal) hypothesis at p=%f\n',p);
        else
            fprintf(' reject H0 (normal) hypothesis at p=%f\n',p);
        end
        break
        
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