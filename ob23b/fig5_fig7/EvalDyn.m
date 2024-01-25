classdef EvalDyn 
    methods     (Static=true)

        function [] = dyn_avg(r_avg, f_avg, sigma_avg, sigma_norm_avg, flag_success, num_gen)
            fprintf('eval.dyn_avg() \n');
            num_success = sum(flag_success); % count successful runs
            N_G = length(r_avg)-1;
            
            r_avg(r_avg==0) = nan; % r_avg initialized with 0 (values added during sim.); "0" no values were written => "nan"
            mask_nan = ones(N_G+1,1); mask_nan(isnan(r_avg)) = nan; % mask_nan sets not written elements as "nan" (N_G >(!) largest occuring generation)
            %title_str = [title_str, ', $P_S$=',num2str(num_success/TRIALS, 2)]; % evaluate P_s
            fprintf('\t min(num_gen)=%i, max(num_gen)=%i \n', min(num_gen), max(num_gen));

            % Average values over num_success
            r_avg = r_avg/num_success .* mask_nan;
            f_avg = f_avg/num_success .* mask_nan;
            sigma_avg = sigma_avg/num_success .* mask_nan;
            sigma_norm_avg = sigma_norm_avg/num_success .* mask_nan;

            %c = linspecer(4)+0.15; c(c>1)=1;
            %fig = figure;
                hold on;
                g=0:1:N_G; 
                    plot(g, r_avg, 'k-', 'DisplayName', '$R$: avg');
                    %plot(g, f_avg, 'color', c(2,:), 'DisplayName', '$f$');
                    %plot(g, sigma_avg, 'color', c(3,:), 'DisplayName', '$\sigma$');
                    %plot(g, sigma_norm_avg, 'color', c(4,:), 'DisplayName', '$\sigma^*$');
                %set(gca, 'YScale', 'log');
                %legend('location', 'best');
                %yticks([1e-6, 1e-4, 1e-2, 1, 1e2, 1e4]);
        end % function
        
        function [] = dyn_norm_avg(r_g_t, flag_success, PLOT_ALL_R, PLOT_NORM_G)
            %% Evaluate normalized data
            fprintf('eval.dyn_norm_avg() \n');
            num_success = sum(flag_success==1);
            r_g_t = r_g_t(:, flag_success==1); % get all succ. runs, discard other
            num_intp = 101; % set interp. point number
            ids = 20; %floor(num_intp/10); % use last (end-ids):end for lin. regression
            data(1:num_success) = struct;

            for i=1:num_success
                id_notNan = ~isnan(r_g_t(:, i));    % get all valid elements (not nan)
                data(i).r = r_g_t(id_notNan,i);     % extract all r(g) valid nums
                data(i).n_data = length(data(i).r); % get number of valid data points
                data(i).g = linspace(0, 1, data(i).n_data)';    % create normalized range for g
                % data(i).m = (log10(data(i).r(end)) - log10(data(i).r(end-ids))) / ids;  % difference quotient?
                c = polyfit((0:1:ids)', log10(data(i).r(end-ids:end)), 1); % get slope in semi-log plot
                data(i).m = c(1); % save data
            end
            m_const = mean([data.m]); % mean of slopes is target slope for avg. curve

            %% Interpolate data
            g_intp = linspace(0,1,num_intp)';   % define g in [0,1] with num_intp
            r_g_t_intp = zeros(num_intp, num_success);
            for i=1:num_success
                r_g_t_intp(:,i) = interp1(data(i).g, data(i).r, g_intp); % => avg. in lin-space
                % r_g_t_intp(:,i) = interp1(data(i).g, log10(data(i).r), g_intp); % => avg. in log-space
            end

            %% Mean curve and slope
            r_g_mean = mean(r_g_t_intp, 2); % average in normalized g (simulation time) space

            % m_intp = (log10(r_g_mean(end)) - log10(r_g_mean(end-ids))) / (g_intp(end) - g_intp(end-ids)); % difference quotient
            c = polyfit(g_intp(end-ids:end), log10(r_g_mean(end-ids:end)), 1); % get slope in normalized space (semilogy)
            m_intp = c(1);
            n_gen_intp = m_intp/m_const; % get number of generations such that mean target slope m_const is reached

            %% Plot interpolated function
            if PLOT_ALL_R == 1
                hold on
            %     for i=1:num_success
            %         % REMARK: data(i).n_data = num of data points; n_data - 1 is the generation value
            %         if PLOT_NORM_G==1
            %             plot(g_intp, r_g_t_intp(:,i), 'b-.'); % plot interpolation curves (reg.)
            %         else
            %             plot(g_intp*(data(i).n_data-1), r_g_t_intp(:,i), 'b-.'); % plot interpolation curves (reg.)     
            %         end
            %     end
                if PLOT_NORM_G==1
                    plot(g_intp, r_g_mean, 'r-', 'DisplayName', '$R$: norm_avg'); % reg g    
                else
                    plot(g_intp*n_gen_intp, r_g_mean, 'r-', 'DisplayName', '$R$: norm_avg') % reg g
                end
            end

            %% Test
            x_test = g_intp(end-ids:end)*n_gen_intp; % regular g space
            y_test = log10(r_g_mean(end-ids:end));   % semilogy
            c = polyfit(x_test, y_test, 1);          % get slope in reg. space
            fprintf('\t -----TEST----- \n');
            fprintf('\t m_const=%d \n', m_const);
            fprintf('\t m_avg=%d \n', c(1));
        end %function
        
        function [] = dyn_avg_glob_attr(r_g_t, flag_success)  
            
            PLOT_ALL_LIN = 0; % plot all linear regressions
            pow_target = -1;  % start lin. reg. at R^pow_target
            pow_target_end = -4; % end lin. reg.
            R_target = 10^pow_target;
            gen = [0:1:size(r_g_t, 1)-1]'; % create gen vector
            
            %% Evaluate average conv. rate of global attr.
            hold on;
            fprintf('eval.dyn_avg_glob_attr() \n');
            num_success = sum(flag_success==1); % get only succ. runs
            r_g_t = r_g_t(:, flag_success==1);  % get only succ. runs
            
            % y=k*gen+d
            k = nan*zeros(1, num_success); % slope of succ. run
            d = nan*zeros(1, num_success); % d of succ. run
            g_endvalue = nan*zeros(1, num_success);
            
            for i=1:num_success
                log_get = r_g_t(:,i)<R_target; % select values below target
                y = log10(r_g_t(log_get, i));  % log10 of R-values as y
                x = gen(log_get);              % gen as x
                c = polyfit(x, y, 1);
                k(i) = c(1); 
                d(i) = c(2);
                if PLOT_ALL_LIN==1
                    px = linspace(x(1)-1, x(end), 101);
                    py = k(i)*px + d(i);
                    plot(px, 10.^py, 'b--'); 
                end
                
                % get average of last generation before convergence => plausibility of linear avg.
                g_endvalue(i) = x(end); 
            end
               
            assert(any(~isnan(k))); assert(any(~isnan(d))); % check
            g_start = (pow_target-d)./k; % get g-value of lin. function (each i) for which target reached (horiz. position)
            k_avg = mean(k);             % avg. slope
            g_avg = mean(g_start);       % avg. crossing point
            d_avg = pow_target - k_avg*g_avg; % translates to avg. d-value, y=kx+d
            g_end_avg = mean(g_endvalue);
            
            %% Plot avg.
            g_endvalue = (pow_target_end-d_avg)/k_avg; % choose some g-end-value
            px = linspace(g_avg, g_endvalue, 101);     % g-range
            py = k_avg*px + d_avg;                % transform back into 10^(.) as plot is taking log(data)        
            plot(px, 10.^py, 'k-');
                %xline(g_end_avg, 'k:', 'LineWidth', 2); % additional info for plausibility of linear avg.
            fprintf('\t R_tgt=%d \n', R_target);
            fprintf('\t k_avg=%d \n', k_avg);
            fprintf('\t g_avg=%d \n', g_avg);
            fprintf('\t d_avg=%d \n', d_avg);
        end%function
        
        function [] = dyn_avg_btw_target(r_g_t, flag_success)  
            
            gen = [0:1:size(r_g_t, 1)-1]'; % create gen vector
            
            hold on;
            fprintf('eval.dyn_avg_btw_target() \n');
            num_success = sum(flag_success==1); % get only succ. runs
            r_g_t = r_g_t(:, flag_success==1);  % get only succ. runs
            log_r_g_t = log10(r_g_t); % giving power of 10
            
            pow_R_max = max(log_r_g_t,[],'all'); % find largest value
            pow_R_min = min(log_r_g_t,[],'all'); % find smallest value
            
            num_intervals = 100;
            pow_range = linspace(pow_R_min, pow_R_max, num_intervals+1)'; % create range of powers btw. min and max value
            pow_counts = zeros(length(pow_range)-1,1);  % count how many elements are added => avg.
            pow_sum = zeros(length(pow_range)-1,1);     % sum of generations between [pow_R(j)+pow_R(j+1)]
            
            for i=1:num_success
                for j=1:num_intervals
                    pow_r_low = pow_range(j);    % get interval
                    pow_r_high = pow_range(j+1); % get interval
                    
                    get = pow_r_low<=log_r_g_t(:,i) & log_r_g_t(:,i)<=pow_r_high;   % get ids, such that value btw. [pow_R(j)+pow_R(j+1)]
                    gen_in_range = gen(get);                                        % get generation values
                    pow_sum(j) = pow_sum(j)+sum(gen_in_range);                      % sum up generation values
                    pow_counts(j) = pow_counts(j)+length(gen_in_range);             % count num. elements
                end      
            end
               
            g_avg = pow_sum ./ pow_counts;                           % average
            pow_avg = 0.5 * (pow_range(1:end-1) + pow_range(2:end)); % values defined for interval, plot point at center value
            
            %% Plot avg.     
            plot(g_avg, 10.^pow_avg, 'g.');
        end%function
    end %methods
end


