function [fig] = fct_plot(pc_input, pc_res, pc_data, init_stop, variation_options, fit, THETA)
    

    % Write Data to text
    %data_cell = ['$f_{st}$=',num2str(pc_res.f, '%.2e'), '; ', '$f_{ev}$=',num2str(pc_res.feval, '%.2e'), '; ', ...
     %   '$\mu_{med}$=',num2str(pc_res.m_median) ];
     %   %     , '; ' ,'$\mu_{q1}$=',num2str(pc_res.m_q1), '; ' ,'$\mu_{q2}$=',num2str(pc_res.m_q2)];
    
    current = variation_options(1);
    num_var = variation_options(2);
    var_type = variation_options(3);
    gen = 1:length(pc_res.r_g);
    N = fit.N;
    
    if current == 1 && (var_type == 1 || var_type == 2)
        fig = figure;
        tiledlayout(num_var,2,'TileSpacing','compact','Padding','compact');
    elseif current == 1 && var_type == 3
        figure;
        tiledlayout(num_var+1,1,'TileSpacing','compact','Padding','compact');    
    end
    
    if var_type ~=3 || (var_type == 3 && current == 1)
        nexttile; hold on;
        %subplot(num_var,2,2*current-1); hold on;
            %plot(gen, pc_res.r_g, 'r-', 'DisplayName', '$R$');
            if ~strcmpi(fit.name, 'Random')
                plot(gen, pc_res.f_g, 'b-.','DisplayName', '$f$');
            end
            myc = [1,0,0]; %[0.39,0.83,0.07];
            plot(gen, pc_res.sigma_g, ':', 'color', myc,'DisplayName', '$\sigma$');
            
            %% Enable to show sigma*
            % if ~strcmpi(fit.name, 'Random')
            %     plot(gen, pc_res.sigma_g./pc_res.r_g*N, 'm-', 'DisplayName', '$\sigma^*$');
            % end
            
            %% Enable to show signzero
            % signzero = sqrt( sqrt(N*(8*e_vartheta_a_b(THETA,1,0).^2.*pc_res.mu_g.^2 + N)) - N);
            % if ~strcmpi(fit.name, 'Random')
            %     plot(gen, signzero, 'k--', 'DisplayName', '$\sigma^*_0$');
            % end

            plot(gen, pc_res.mu_g, 'k', 'DisplayName', '$\mu$')

            %xlabel('$g$'); 
            ylabel([fit.label, ' ($N$=', num2str(N), ')']);
            set(gca, 'YScale', 'log');
            %xl = xlim;
            %yticks(init_stop.yticks);
            %ylim(init_stop.ylim);

            %if current == 1, legend('location','best','NumColumns',2); end

            yticks([10.^(-20:4:14)]);
    end
        
        %% Choose data for plot
    %     subplot(num_var,2,2*current); hold on;
        nexttile; hold on;
        switch pc_input.pc_method
            
            case 'ROSAES'
                % --- dual axes ---
                %negid = (pc_data.dF<0);
                %posid = ~negid;
                %yyaxis left
                %plot(gen(negid),pc_data.dF(negid)); xlim(xl);
                %set(gca,'Yscale','log')
                %yyaxis right
                %plot(gen(posid),pc_data.dF(posid)); xlim(xl);
                
                plot(gen,pc_data.dF); 
                xl = xlim; xlim([0,xl(2)]);
                y_plot = max(pc_data.dF);
                if y_plot>0
                    ylim([-4*y_plot, 1.2*y_plot])
                end
                yline(pc_input.robust_thresh, 'k')
                ylabel('$\overline{\Delta f}$');
                
            case 'PCCMSA'
                pc_data.PH = pc_data.PH(gen);
                plot(gen,pc_data.PH, 'k-'); 
                xl = xlim; xlim([0,xl(2)]);
                ylim([0,1])
                yline(pc_input.alpha_sign, 'm--')
                ylabel('$P_H$');
                % yticks([0:0.2:1]);  set(gca,'YScale','lin');
                yticks([0.01,0.1,1]); set(gca,'YScale','log'); ylim([0.01,1])
                %ids = find(H_g_t==0);
                %gen = gen;
                %plot(gen(ids), r_g_t(ids), '.'); 
                
            case 'APOP'
                plot(gen,pc_data.N_up, 'k-');  
                xl = xlim; xlim([0,xl(2)]);
                ylim([0,1])
                yline(pc_input.apop_thresh, 'm--')
                yticks([0:0.2:1]);
                ylabel('$P_{f}$'); 
                
            case 'INFNOISE'
                plot(gen,pc_data.D2, 'b-'); 
                plot(gen,pc_data.P_INFNOISE, 'k-'); 
                xl = xlim; xlim([0,xl(2)]);
                ylabel('$P_\mathrm{INFNOISE} | D_x^2$');
                if pc_input.pc_on == 1
                    yline(pc_input.alpha_infnoise, 'c-')          
                    %yline(pc_input.ev_infnoise(pc_input.mu_max), '--')
                    %yline(pc_input.ev_infnoise(pc_input.mu_max)-2*pc_input.stddev_infnoise(pc_input.mu_max,N), '--')
                else
                    %yline(pc_input.ev_infnoise(pc_input.mu_min), '--')
                    %yline(pc_input.ev_infnoise(pc_input.mu_min)-2*pc_input.stddev_infnoise(pc_input.mu_min,N), '--')
                end
                    % if pc_input.pc_on == 0
                %     D2_randSel = (1-1/MU);
                %     EV_m_linFct_largeN = (1-1/MU)*(1+1/N*(e_mu_lam_a_b_v2(MU,LAM,1,1)-e_mu_lam_a_b_v2(MU,LAM,2,0)));
                %     yline(D2_randSel)
                %     D2_var_analytic = 2/(MU*N)*(1-1/MU);
                %     yline(D2_randSel+3*sqrt(D2_var_analytic), '--')
                %     yline(D2_randSel-3*sqrt(D2_var_analytic), '--')
                %     yline(EV_m_linFct_largeN, 'r--')
                %     D2_var_meas = std(pc_data.D2).^2;
                %     ratio = D2_var_meas/D2_var_analytic
                % end
            case 'QGAIN'
                mean_Q = mean(pc_data.Q, 'omitnan');
                std_Q = std(pc_data.Q, 0, 'omitnan');
                
                plot(gen,pc_data.Q);  
                xl = xlim; xlim([0,xl(2)]);
                if strcmpi(name, 'Sphere')
                    sigma_norm = mean(sigma_g./r_g*N, 'omitnan');
                    noise = 0;
                    Q_sphere = sphere_qgain(MU, LAM, N, sigma_norm, noise);
                    yline(mean_Q, 'k--');
                    yline(-Q_sphere, 'r-');
                end
                %ylim([0,1])
                ylabel('$\Delta f^*$');  
                
            case 'PSA'
                plot(pc_data.norm2_p_dm, 'b-'); 
                plot(pc_data.norm2_p_dS, 'r-'); 
                plot(pc_data.norm2_p_dm+pc_data.norm2_p_dS, 'k-');  
                yline(pc_input.alpha_psa, 'm--');
                ylabel('$||\mathbf{p}_\theta||^2$'); 

            case 'LENGTHZ'
                plot(vecnorm(pc_data.Z,2,1).^2, 'k-');  
                yline(1);
                yline(1-2*sqrt(2/fit.N));
                yline(1+2*sqrt(2/fit.N));
                ylabel('$z^2$'); 

            case 'RANK'
                pc_data.P_RANK = pc_data.P_RANK(gen);
                plot(gen,pc_data.P_RANK, '-'); 
                yline(pc_input.alpha_rank, 'k--')
                xl = xlim; xlim([0,xl(2)]);
                ylim([0,1])
                yline(pc_input.alpha_rank)
                ylabel('$P_\mathrm{RANK}$');
                yticks([0:0.2:1]);
        end   
        if pc_input.annotate==1
            x0 = xlim;
            y0 = ylim;
            text(x0(1)+10,y0(2)*0.95,data_cell)
        end
        if pc_input.annotate==1
            data_cell = fieldnames((pc_input))+string(": "+struct2cell(pc_input));
            annotation('textbox',[0.85,0.68,0.3,0.3],'String',data_cell,'Interpreter','None','FontSize',8,'LineStyle','None')
        end
        myfigstyle(gcf,14,num_var*3,8,8);
end