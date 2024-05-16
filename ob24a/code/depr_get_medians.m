function [] = get_medians(r_g_t, sigma_g_t, flag_success, vec_select, N)

    get_ids = [flag_success,~flag_success];

%     phi_g = r_g(1:end-1).^2-r_g(2:end).^2;
%     phi_g = [phi_g;nan];
%     phi_norm_g = phi_g*N./(2*r_g.^2);
    
    %% Plot all dynamics
    figure; hold on;
    %subplot(1,3,1); hold on;
        %legend('AutoUpdate','off');
        [N_G, TRIALS] = size(r_g_t);
        for i=1:TRIALS
            if flag_success(i)==1, c = [0.6,0.6,0.6]; else, c = [0.9290, 0.6940, 0.1250]; end
            plot(0:1:N_G-1, r_g_t(:,i), 'color', c);
        end
        %legend('AutoUpdate','on');
        
    for i = 1:2
    vec_select = flag_success;
   
    %% Get selected runs
    r_g = median(r_g_t(:,logical(vec_select)), 2, 'includenan');
    sigma_g = median(sigma_g_t(:,logical(vec_select)), 2, 'includenan');
    sigma_norm_g = median(sigma_g_t(:,logical(vec_select))*N./r_g, 2, 'includenan');

    %% Remove nan
    r_g = r_g(~isnan(r_g));
    sigma_g = sigma_g(~isnan(sigma_g));
    sigma_norm_g = sigma_norm_g(~isnan(sigma_g));
    gen = 0:1:length(r_g)'-1;       
        
    %% Plot median
    plot(gen, r_g, 'k-', 'DisplayName', 'median($R$)');
    xlabel('$g$'), ylabel('$R(g)$');
    set(gca, 'YScale', 'log');
    myfigstyle(gcf, 16, 10, 10, 10);
    
    %% Sigma*, phi*
    figure; hold on;
    %subplot(1,3,2); hold on;
        for i=1:TRIALS
            if flag_success(i)==1, c = [0.6,0.6,0.6]; else, c = [0.9290, 0.6940, 0.1250]; end
            plot(0:1:N_G-1, sigma_g_t(:,i)./r_g_t(:,i)*N, 'color', c);
        end
        plot(gen, sigma_norm_g, 'k-', 'DisplayName', 'median($\sigma^*$)');
        xlabel('$g$'); ylabel('$\sigma^*$-dynamics');
        set(gca, 'YScale', 'log');
%     subplot(1,3,3); hold on;
%         plot(gen, phi_norm_g, 'k-');
%         xlabel('$\sigma^*$'); ylabel('$\varphi^*_R$-dynamics');
    figure; hold on;
    %subplot(1,3,3); hold on;
        plot(sigma_norm_g, r_g, 'k-');
        xlabel('$\sigma^*$'); ylabel('$R(\sigma^*)$-dynamics');
        set(gca, 'YScale', 'log');
    myfigstyle(gcf, 32, 10, 10, 10);
    
    %% DEBUG SIGMA_NORM
    %figure; plot(sigma_g_t*N./r_g_t)
end