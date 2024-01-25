function [] = get_medians(r_g_t, sigma_g_t, flag_success, vec_select, N)

    %assert(all(size(vec_select)==size(r_g_t(1,:))));
    vec_select = reshape(vec_select, 1, length(vec_select));

    %% Get selected runs
    r_g = median(r_g_t(:,logical(vec_select)), 2, 'includenan');
    sigma_g = median(sigma_g_t(:,logical(vec_select)), 2, 'includenan');
    sigma_norm_g = median(sigma_g_t(:,logical(vec_select))*N./r_g, 2, 'includenan');
    gen = 0:1:length(r_g)'-1;
    
    %% Plot all dynamics

    figure;  hold on;
    cloc = colormap(cool(10));
    cconv = colormap(summer(10));
        %legend('AutoUpdate','off');
        [N_G, TRIALS] = size(r_g_t);
        for i=1:TRIALS
            if flag_success(i)==1
                c = cconv(randi(10),:); %[0.6,0.6,0.6]; 
            else
                c = cloc(randi(10),:); %[0.9290, 0.6940, 0.1250]; 
            end
            plot(0:1:N_G-1, r_g_t(:,i), 'color', c);
        end
        %legend('AutoUpdate','on');
        %% Plot median
        plot(gen, r_g, 'k-', 'DisplayName', 'median($R$)');
        xlabel('$g$'), ylabel('$R(g)$');
        set(gca, 'YScale', 'log');
    myfigstyle(gcf, 10, 7, 9, 9);
        
    figure; hold on;
        for i=1:TRIALS
            if flag_success(i)==1
                c = cconv(randi(10),:); %[0.6,0.6,0.6]; 
            else
                c = cloc(randi(10),:); %[0.9290, 0.6940, 0.1250]; 
            end
            plot(0:1:N_G-1, sigma_g_t(:,i)./r_g_t(:,i)*N, 'color', c); %
        end
        plot(gen, sigma_norm_g, 'k-');
        xlabel('$g$'); ylabel('$\sigma^*(g)$');
    myfigstyle(gcf, 10, 7, 9, 9);
    
    % figure; 
    % plot(sigma_norm_g, r_g, 'k-');
    % xlabel('$\sigma^*(g)$'); ylabel('$R(g)$');
    % set(gca, 'YScale', 'log');
    % myfigstyle(gcf, 8, 5, 9, 9);
    
end