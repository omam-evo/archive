function [r_g,sigma_g,sigma_norm_g] = get_means(r_g_t, sigma_g_t, flag_success, vec_select, N, OPTIONS)

    %assert(all(size(vec_select)==size(r_g_t(1,:))));
    vec_select = reshape(vec_select, 1, length(vec_select));

    %% Get selected runs
    r_g = mean(r_g_t(:,logical(vec_select)), 2, 'includenan');
    sigma_g = mean(sigma_g_t(:,logical(vec_select)), 2, 'includenan');
    sigma_norm_g = mean(sigma_g_t(:,logical(vec_select))*N./r_g, 2, 'includenan');
    gen = (1:1:length(r_g))';

    cloc = repmat([0.8,0.8,0.8],10,1); %colormap(cool(10));
    cconv = repmat([0.6,0.6,0.6],10,1); %colormap(summer(10));
    
    %% Plot all dynamics
    if any(strcmpi(OPTIONS, 'dyn'))
        figure;  hold on;
            %legend('AutoUpdate','off');
            [N_G, TRIALS] = size(r_g_t);
            if any(strcmpi(OPTIONS, 'all'))
                for i=1:TRIALS
                    if flag_success(i)==1
                        c = cconv(randi(10),:); %[0.6,0.6,0.6]; 
                    else
                        c = cloc(randi(10),:); %[0.9290, 0.6940, 0.1250]; 
                    end
                    plot(gen, r_g_t(:,i), 'color', c);
                end
            end
            %legend('AutoUpdate','on');
            %% Plot mean
            plot(gen, r_g, 'k-', 'DisplayName', 'mean($R$)');
            xlabel('$g$'), ylabel('$R(g)$');
            set(gca, 'YScale', 'log');
        myfigstyle(gcf, 10, 7, 9, 9);
    end
        
    if any(strcmpi(OPTIONS, 'sign_g'))
        figure; hold on;
            if any(strcmpi(OPTIONS, 'all')) 
                for i=1:TRIALS
                    if flag_success(i)==1
                        c = cconv(randi(10),:); %[0.6,0.6,0.6]; 
                    else
                        c = cloc(randi(10),:); %[0.9290, 0.6940, 0.1250]; 
                    end
                    plot(0:1:N_G-1, sigma_g_t(:,i)./r_g_t(:,i)*N, 'color', c); %
                end
            end
            plot(gen, sigma_norm_g, 'k-');
            xlabel('$g$'); ylabel('$\sigma^*(g)$');
        myfigstyle(gcf, 10, 7, 9, 9);
    end
    if any(strcmpi(OPTIONS, 'sig_g'))
        figure; hold on;
            if any(strcmpi(OPTIONS, 'all')) 
                for i=1:TRIALS
                    if flag_success(i)==1
                        c = cconv(randi(10),:); %[0.6,0.6,0.6]; 
                    else
                        c = cloc(randi(10),:); %[0.9290, 0.6940, 0.1250]; 
                    end
                    plot(0:1:N_G-1, sigma_g_t(:,i), 'color', c); %
                end
            end
            plot(gen, sigma_g, 'k-');
            xlabel('$g$'); ylabel('$\sigma(g)$');
        myfigstyle(gcf, 10, 7, 9, 9);
    end
    if any(strcmpi(OPTIONS, 'sign_r'))
        figure; 
        plot(sigma_norm_g, r_g, 'k-');
        xlabel('$\sigma^*(g)$'); ylabel('$R(g)$');
        set(gca, 'YScale', 'log');
        myfigstyle(gcf, 8, 5, 9, 9);
    end
    
end