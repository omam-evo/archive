%% Infos
% check plausibility on random function:  sum(pc_data.PH(:,1)<=pc_input.alpha_sign)/G

function [PH_g, H_g, T] = pc_pccmsa(pc_input, g, PH_g, H_g, F_g, DEBUG_G)
        g0 = g-(pc_input.L-1);

        %% Get values for lin. regression
        xi = (g0:g)'; %AR_x(g-g0:g);
        yi = F_g(g0:g);
        %% Log-values for lin. regression
%         f_min = min(AR_f(g-g0:g));
%         if f_min<=0
%            yi = log10(AR_f(g-g0:g)-f_min+1e-6);
%         else
%            yi = log10(AR_f(g-g0:g));
%         end         
        %% Linear: y = a*x+b
        %p = polyfit(xi,yi,1); 
        xm = mean(xi);
        ym = mean(yi);
        a = sum((xi-xm).*(yi-ym))/sum((xi-xm).^2); % measured quantity => to be investigated
        b = ym-a*xm; %measured

        std_error_a = sqrt(sum((yi-ym).^2)/(pc_input.L-2)/sum( (xi-xm).^2 ));

        %% Probability of a random slope "a" to be below measured
        a_H0 = 0;
        T_a = (a-a_H0) / std_error_a;
        PH_g(g+1) = tcdf(T_a, pc_input.L-2);
        % PH_g(g+1:g+ids) = nan; %use if future L-values are nan due to wait (if function not called again)

        % t-value to reach significance level
        t_alpha = tinv(pc_input.alpha_sign, pc_input.L-2); %tcdf(t_alpha,L-2) = alpha_sign
        H_g(g+1) = (a<std_error_a*t_alpha); % compare
        % H_g(g+1:g+ids) = nan;  %use if future L-values are nan due to wait
        

        %% Plot t-Student distribution with L-2 dof
        if ~isnan(DEBUG_G) && g==DEBUG_G
            a_H0 = 0;
            x = linspace(-6*std_error_a, 6*std_error_a, 1001);
            t = (x-a_H0)/std_error_a;
            figure; 
                subplot(1,2,1); hold on; 
                    plot(xi,yi);
                    plot(xi, a*xi+b, '-');
                    xlabel('$g$'); ylabel('$f(g)$');
                    title(['$\hat{a}$=',num2str(a), ', $s_{\hat{a}}$=',num2str(std_error_a), ', $\hat{b}$=',num2str(b)])
                subplot(1,2,2); hold on; legend;
                    plot(t, tpdf(t,pc_input.L-2), 'k-', 'DisplayName','$t$-PDF (L-2)');
                    xlabel('$t$'); ylabel('PDF');
                    xline(t_alpha, 'b--', 'DisplayName', ['$P_H=',num2str(pc_input.alpha_sign),'$'])
                    xline(T_a, 'r', 'DisplayName', ['$P_H=',num2str(PH_g(g+1)),'$']) 
                    %title(['$T_a$=$(a-a_{H_0}) / s_a$=',num2str(T_a)])
                %sgtitle(['H_g(g=',num2str(g),')=',num2str(H_g(g))]);
            myfigstyle(gcf,30,10,11,11);
            pause
        end
 
        x = PH_g(g+1);
        % if strcmpi(pc_input.pid_method, 'pid_v1') 
            if x < pc_input.alpha_sign
                T = 1;
            elseif x == pc_input.alpha_sign
                T = 0;
            else
                T = -1;
            end
        % end
        % elseif strcmpi(pc_input.pid_method, 'pid_v2')
            % if x < pc_input.alpha_sign
            %     k = -1/pc_input.alpha_sign;
            %     d = 1;
            %     T = k*x+d;
            % elseif x == pc_input.alpha_sign
            %     T = 0;
            % else
            %     k = -1/(1-pc_input.alpha_sign);
            %     d = pc_input.alpha_sign/(1-pc_input.alpha_sign);
            %     T = k*x+d;
            % end        
        % end
        
end

