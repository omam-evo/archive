classdef EvalPhi
    methods     (Static=true)

        function [phi_2_num, phi_2_B1, phi_2_B2, phi_2_L1, phi_2_L2] = get_phi2_i(ras_tools, MU, LAM, A, ALPHA, SIGMA_LIST, Y_MAT, phi_2_mean, phi_2_stderr, FLAG_SIGMA_NORM, COMP_ID, SAVEPATH, OPTIONS)
            
            N = size(Y_MAT,1);
            NUM_SIG = length(SIGMA_LIST);
            NUM_POS = size(Y_MAT,2);

            %% INIT PHI
            Phi = Phi_v2;

            %phi_1_A3 = zeros(N, NUM_SIG, NUM_POS);
            %phi_1_C_k = zeros(N, NUM_SIG, NUM_POS);
            %phi_1_C_f = zeros(N, NUM_SIG, NUM_POS);

            %% INIT PHI_2
            phi_2_num = zeros(NUM_SIG, NUM_POS);
            e_10 = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
            e_11 = e_mu_lam_a_b_v2(MU, LAM, 1, 1);       % only needed for phi_2
            e_20 = e_mu_lam_a_b_v2(MU, LAM, 2, 0);       % only needed for phi_2
            phi_2_B1 = zeros(N, NUM_SIG, NUM_POS);
            phi_2_B2 = zeros(N, NUM_SIG, NUM_POS);
            phi_2_L1 = zeros(N, NUM_SIG, NUM_POS);
            phi_2_L2 = zeros(N, NUM_SIG, NUM_POS);
            
            if any(strcmpi(OPTIONS, 'L1_TERMS'))
                phi_2_L1_T1 = nan*phi_2_L1;
                phi_2_L1_T2 = nan*phi_2_L1;
                phi_2_L1_T3 = nan*phi_2_L1;
            end

            %% RUN APPROXIMATION
            for p=1:NUM_POS
                Y_0 = Y_MAT(:, p);
                Y_i = Y_0(COMP_ID); Y_jni = Y_0(COMP_ID ~= 1:N);

                for s=1:NUM_SIG

                    if FLAG_SIGMA_NORM==0, sigma = SIGMA_LIST(s); else, sigma = SIGMA_LIST(s)*norm(Y_0)/N; end

                    %% Phi 1
                    %[phi, ~, ~] = Phi.get_phi_1_viaC_Di(ras_tools, C_MU_LAM, A, ALPHA, sigma, Y_0, 0);
                    %phi_1_C_k(:, s, p) = phi;
                    %[phi, ~, ~] = Phi.get_phi_1_viaC_Di(ras_tools, C_MU_LAM, A, ALPHA, sigma, Y_0, 1);
                    %phi_1_C_f(:, s, p) = phi;
                    %[phi, ~, ~] = Phi.get_phi_1(ras_tools, MU, LAM, A, ALPHA, sigma, Y_0, 1);
                    %phi_1_A3(:, s, p) = phi;

                    %% Phi 2
                    phi_2_B1(:, s, p) = Phi.get_phi_2(ras_tools, MU, LAM, A, ALPHA, sigma, Y_0, e_11, e_20);
                    phi_2_B2(:, s, p) = Phi.get_phi_2(ras_tools, MU, LAM, A, ALPHA, sigma, Y_0, 0, 0);
                    
                    if any(strcmpi(OPTIONS, 'L1'))
                        phi_2_L1(:, s, p) = Phi.get_phi_2_lgPop(ras_tools, MU, LAM, A, ALPHA, sigma, Y_0, {'L1'});
                    end
                    if any(strcmpi(OPTIONS, 'L1_TERMS'))
                        [phi_2_L1(:, s, p),phi_2_L1_T1(:, s, p),phi_2_L1_T2(:, s, p),phi_2_L1_T3(:, s, p)] = Phi.get_phi_2_lgPop(ras_tools, MU, LAM, A, ALPHA, sigma, Y_0, {'L1'});
                    end
                    if any(strcmpi(OPTIONS, 'L2'))
                        phi_2_L2(:, s, p) = Phi.get_phi_2_lgPop(ras_tools, MU, LAM, A, ALPHA, sigma, Y_0, {'L2'});
                    end
                    if any(strcmpi(OPTIONS, 'NUM')) 
                        phi_2_num(s, p) = Phi.get_phi_2_num_lgPop(ras_tools, A, ALPHA, MU, LAM, sigma, Y_0, Y_i, Y_jni);
                    end
                end
            end

            %% Phi Plots
            if any(strcmpi(OPTIONS, 'PLOT'))
                if NUM_POS==1; num_row = 1; else; num_row = 2; end
                num_col = ceil(NUM_POS/2);
                figure; tiledlayout(num_row, num_col, 'Padding', 'none', 'TileSpacing', 'compact'); 
                for p=1:NUM_POS
                    Y_0 = Y_MAT(:, p);
                    Y_i = Y_0(COMP_ID); Y_jni = Y_0(COMP_ID ~= 1:N);

                    nexttile; hold on;
                    if p==1, legend('location', 'southwest'); end

                    %% Title and labels
                    %Y_str = []; %,', Y_j=[', num2str(Y_inj(1),3),', ',num2str(Y_inj(2),3),', ...,',num2str(Y_inj(end),3), ']']; 
                    title_str = ['(',num2str(MU),',',num2str(LAM),')-ES','~$N\!=\!',num2str(N),'~A\!=\!',num2str(A),'~R\!=\!',num2str(norm(Y_0)),'~y_i\!=\!', num2str(Y_i,3), '~(i\!=\!',num2str(COMP_ID),')$'];
                    if FLAG_SIGMA_NORM==0; str_sig_label = '$\sigma$'; else, str_sig_label = '$\sigma^*$'; end
                    title(title_str);
                    xlabel(str_sig_label); ylabel('$\varphi_i^{\mathrm{II}}$');

                        %% Phi experimental 
                        if any(strcmpi(OPTIONS, 'SIM'))
                            errorbar(SIGMA_LIST, phi_2_mean(COMP_ID, :, p), phi_2_stderr(COMP_ID, :, p), 'k.', 'DisplayName', 'SIM');
                        end
                        if any(strcmpi(OPTIONS, 'NUM')) % show numeric  
                            plot(SIGMA_LIST, phi_2_num(:, p), 'k:', 'DisplayName', 'NUM LP');
                        end
                        %% Phi Approx
                        if any(strcmpi(OPTIONS, 'B1'))
                            plot(SIGMA_LIST, phi_2_B1(COMP_ID, :, p), 'r:', 'DisplayName', 'B1');
                        end
                        if any(strcmpi(OPTIONS, 'B2'))
                            plot(SIGMA_LIST, phi_2_B2(COMP_ID, :, p), 'b--', 'DisplayName', 'B2');
                        end
                        if any(strcmpi(OPTIONS, 'L1'))
                            plot(SIGMA_LIST, phi_2_L1(COMP_ID, :, p), 'g-.', 'DisplayName', 'L1');
                        end
                        if any(strcmpi(OPTIONS, 'L2'))
                            plot(SIGMA_LIST, phi_2_L2(COMP_ID, :, p), 'm--', 'DisplayName', 'L2');
                        end  

                end %NUM_POS
                myfigstyle(gcf, 12*num_col, 6*num_row, 10, 10);
                if any(strcmpi(OPTIONS, 'SAVE'))
                    saveas(gcf,fullfile(SAVEPATH,['phi2','_comp',num2str(COMP_ID),'.fig']));
                end
                
                if any(strcmpi(OPTIONS, 'L1_TERMS'))
                    figure; tiledlayout(num_row, num_col, 'Padding', 'none', 'TileSpacing', 'compact'); 
                    for p=1:NUM_POS  
                        Y_0 = Y_MAT(:, p);
                        Y_i = Y_0(COMP_ID); Y_jni = Y_0(COMP_ID ~= 1:N);

                        nexttile; hold on;
                        if p==1, legend('location', 'southwest'); end

                        %% Title and labels
                        title_str = ['(',num2str(MU),',',num2str(LAM),')-ES','~$N\!=\!',num2str(N),'~A\!=\!',num2str(A),'~R\!=\!',num2str(norm(Y_0)),'~y_i\!=\!', num2str(Y_i,3), '~(i\!=\!',num2str(COMP_ID),')$'];
                        if FLAG_SIGMA_NORM==0; str_sig_label = '$\sigma$'; else, str_sig_label = '$\sigma^*$'; end
                        title(title_str);
                        xlabel(str_sig_label); ylabel('Terms of $\varphi_i^{\mathrm{II}}$');
                            yline(1,'k--', 'DisplayName', 'Value: 1');
                            plot(SIGMA_LIST, phi_2_L1_T1(COMP_ID, :, p), '--', 'DisplayName', 'Term T1');
                            plot(SIGMA_LIST, phi_2_L1_T2(COMP_ID, :, p), '--', 'DisplayName', 'Term T2');
                            plot(SIGMA_LIST, phi_2_L1_T3(COMP_ID, :, p), '--', 'DisplayName', 'Term T3');

                    end
                    myfigstyle(gcf, 12*num_col, 6*num_row, 10, 10);
                end
                if any(strcmpi(OPTIONS, 'SAVE'))
                    saveas(gcf,fullfile(SAVEPATH,['phi2','_terms_comp',num2str(COMP_ID),'.fig']));
                end
            end % DO PLOT

        end %function get_phi_2_plot
        
        function [phi_R2_mean,phi_R2_B1,phi_R2_B2,phi_R2_L1,phi_R2_L2] = get_phi2_R2(ras_tools, MU, LAM, A, ALPHA, SIGMA_LIST, Y_MAT, phi_R2_mean, phi_R2_stderr, phi_2_B1, phi_2_B2, phi_2_L1, phi_2_L2, FLAG_NORM_PHI, FLAG_SIGMA_NORM, SAVEPATH, OPTIONS)
            
            N = size(Y_MAT,1);
            NUM_SIG = length(SIGMA_LIST);
            NUM_POS = size(Y_MAT,2);
            
            phi_R2_B1 = reshape(sum(phi_2_B1,1),NUM_SIG,NUM_POS);
            phi_R2_B2 = reshape(sum(phi_2_B2,1),NUM_SIG,NUM_POS);
            phi_R2_L1 = reshape(sum(phi_2_L1,1),NUM_SIG,NUM_POS);
            phi_R2_L2 = reshape(sum(phi_2_L2,1),NUM_SIG,NUM_POS);
            
            phi_R_SPH = nan*zeros(NUM_SIG,NUM_POS);
            
            if NUM_POS==1; num_row = 1; else; num_row = 2; end
            num_col = ceil(NUM_POS/2);
            figure; tiledlayout(num_row, num_col, 'Padding', 'none', 'TileSpacing', 'compact'); 
            for i=1:NUM_POS
                
                %% R
                R = norm(Y_MAT(:,i));  
                
                %% NORM
                if FLAG_NORM_PHI == 1
                    R = norm(Y_MAT(:,i));
                    phi_R2_reg2norm = N/(2*R^2);
                    phi_R2_mean(:,i) = phi_R2_mean(:,i)*phi_R2_reg2norm;
                    phi_R2_stderr(:,i) = phi_R2_stderr(:,i)*phi_R2_reg2norm;
                    phi_R2_B1(:,i) = phi_R2_B1(:,i)*phi_R2_reg2norm;
                    phi_R2_B2(:,i) = phi_R2_B2(:,i)*phi_R2_reg2norm;
                    phi_R2_L1(:,i) = phi_R2_L1(:,i)*phi_R2_reg2norm;
                    phi_R2_L2(:,i) = phi_R2_L2(:,i)*phi_R2_reg2norm;
                    
                    p = Phi_v2;
                    phi_R_SPH(:,i) = p.get_phi_1_sphere(MU, LAM, N, SIGMA_LIST);
                end
                
                %% Title and labels 
                if any(strcmpi(OPTIONS, 'PLOT'))
                    nexttile; hold on;
                    if i==1, legend('location', 'southwest'); end
                    title_str = ['(',num2str(MU),',',num2str(LAM),')-ES','~$N\!=\!',num2str(N),'~A\!=\!',num2str(A),'~R\!=\!',num2str(R),'$'];
                    title(title_str);
                    if FLAG_SIGMA_NORM==0; xlabel('$\sigma$'), else, xlabel('$\sigma^*$'); end
                    if FLAG_NORM_PHI==1 
                        ylabel('$\varphi_R^{\mathrm{II},*}$');
                    else
                        ylabel('$\varphi_R^{\mathrm{II}}$');
                    end

                    %% Phi experimental
                    if any(strcmpi(OPTIONS, 'SIM'))
                        errorbar(SIGMA_LIST, phi_R2_mean(:,i), phi_R2_stderr(:,i), 'k.', 'DisplayName', 'SIM');
                    end
                    
                    %% Phi Approx
                    if any(strcmpi(OPTIONS, 'B1'))
                        plot(SIGMA_LIST, phi_R2_B1(:,i), 'r-.', 'DisplayName', 'B1');
                    end
                    if any(strcmpi(OPTIONS, 'B2'))
                        plot(SIGMA_LIST, phi_R2_B2(:,i), 'b--', 'DisplayName', 'B2');
                    end
                    if any(strcmpi(OPTIONS, 'L1'))
                        plot(SIGMA_LIST, phi_R2_L1(:,i), 'g-.', 'DisplayName', 'L1');
                    end
                    if any(strcmpi(OPTIONS, 'L2'))
                        plot(SIGMA_LIST, phi_R2_L2(:,i), 'm--', 'DisplayName', 'L2');
                    end
                    
                    %% Sphere
                    if any(strcmpi(OPTIONS, 'SPH'))
                        plot(SIGMA_LIST, phi_R_SPH(:,i), ':', 'color', [0.5,0.5,0.5], 'DisplayName', 'Sphere $\varphi^*$');
                    end
                end
                
            end %NUM_POS
            if any(strcmpi(OPTIONS, 'PLOT'))
                myfigstyle(gcf, 12*num_col, 6*num_row, 10, 10);
            end
            
            if any(strcmpi(OPTIONS, 'SAVE'))
                saveas(gcf,fullfile(SAVEPATH,['phi2','_R2','.fig']));
            end
        end
        
        function [phi_1_num, phi_1_A3, phi_1_C_f] = get_phi1_i_plot(ras_tools, MU, LAM, A, ALPHA, SIGMA_LIST, Y_MAT, phi_1_mean, phi_1_stderr, FLAG_SIGMA_NORM, COMP_ID, SAVEPATH, OPTIONS)
        
        if any(strcmpi(OPTIONS, 'ki'))
            bool_ki = 1; bool_di = 0;
        elseif any(strcmpi(OPTIONS, 'di'))
            bool_di = 1; bool_ki = 0;
        else
            bool_ki = 1; bool_di = 1; %no option given, activate both
        end
            
        N = size(Y_MAT,1);
        NUM_SIG = length(SIGMA_LIST);
        NUM_POS = size(Y_MAT,2);
        
        %% INIT PHI
        Phi = Phi_v2;
        %C_MU_LAM = e_mu_lam_a_b_v2(MU, LAM, 1, 0); % only needed for phi_1 via c_mu_lam
        phi_1_A3 = zeros(N, NUM_SIG, NUM_POS);
        %phi_1_C_k = zeros(N, NUM_SIG, NUM_POS);
        phi_1_C_f = zeros(N, NUM_SIG, NUM_POS);
        phi_1_num = zeros(NUM_SIG, NUM_POS);
                
        %% RUN APPROXIMATION
        for p=1:NUM_POS
            Y_0 = Y_MAT(:, p);
            Y_i = Y_0(COMP_ID); Y_jni = Y_0(COMP_ID ~= 1:N);
        
            for s=1:NUM_SIG
        
                if FLAG_SIGMA_NORM==0, sigma = SIGMA_LIST(s); else, sigma = SIGMA_LIST(s)*norm(Y_0)/N; end
                
                %% Phi 1
                %[phi, ~, ~] = Phi.get_phi_1_viaC_Di(ras_tools, C_MU_LAM, A, ALPHA, sigma, Y_0, 0);
                %phi_1_C_k(:, s, p) = phi;
                
                %% LINEARIZED Q GAIN SOLUTION (METHOD Arnold)
                EVAL_di_term = 1;
                [phi, ~, ~] = Phi.get_phi_1_viaC_Di(ras_tools, e_mu_lam_a_b_v2(MU,LAM,1,0), A, ALPHA, sigma, Y_0, EVAL_di_term);
                phi_1_C_f(:, s, p) = phi;
                %% VARTHETA SOLUTION
                [phi, ~, ~] = Phi.get_phi_1(ras_tools, MU, LAM, A, ALPHA, sigma, Y_0, bool_ki, bool_di);
                phi_1_A3(:, s, p) = phi;
                
                if any(strcmpi(OPTIONS, 'NUM')) % show numeric  
                    Y_i = Y_0(COMP_ID); Y_jni = Y_0(COMP_ID ~= 1:N);
                    phi_1_num(s,p) = Phi.get_phi_1_num(ras_tools, A, ALPHA, MU, LAM, sigma, Y_0, Y_i, Y_jni);
                end                
            end
        end
        
        %% Phi Plots
        if NUM_POS==1; num_row = 1; else; num_row = 2; end
        num_col = ceil(NUM_POS/2);
        figure;
        if NUM_POS > 1
            tiledlayout(num_row, num_col, 'Padding', 'none', 'TileSpacing', 'compact'); 
        end
        for p=1:NUM_POS
            Y_0 = Y_MAT(:, p);
            Y_i = Y_0(COMP_ID); Y_jni = Y_0(COMP_ID ~= 1:N);
  
            if NUM_POS > 1, nexttile; end
            hold on;
            if p==1, legend('location', 'southwest'); end

            %% Title and labels
            %Y_str = []; %,', Y_j=[', num2str(Y_inj(1),3),', ',num2str(Y_inj(2),3),', ...,',num2str(Y_inj(end),3), ']']; 
            title_str = ['(',num2str(MU),',',num2str(LAM),')-ES','~$N\!=\!',num2str(N),'~A\!=\!',num2str(A),'~R\!=\!',num2str(norm(Y_0)),'~y_i\!=\!', num2str(Y_i,3), '~(i\!=\!',num2str(COMP_ID),')$'];
            if FLAG_SIGMA_NORM==0; str_sig_label = '$\sigma$'; else, str_sig_label = '$\sigma^*$'; end
            title(title_str);
            xlabel(str_sig_label); ylabel('$\varphi_i$');
 
                if any(strcmpi(OPTIONS, 'SIM'))
                    errorbar(SIGMA_LIST, phi_1_mean(COMP_ID, :, p), phi_1_stderr(COMP_ID, :, p), 'k.', 'DisplayName', 'Sim.');
                end
                if any(strcmpi(OPTIONS, 'NUM')) % show numeric  
                    plot(SIGMA_LIST, phi_1_num(:, p), 'k:', 'DisplayName', 'Num.');
                end
                if any(strcmpi(OPTIONS, 'A3'))
                    if bool_ki==1 && bool_di==1
                        c = 'r-.';
                    elseif bool_ki==1 
                        c = 'c--';
                    elseif bool_di==1
                        c = 'g:';
                    end
                    plot(SIGMA_LIST, phi_1_A3(COMP_ID, :, p), c, 'DisplayName', '$\varphi_i$');% = 2y_i\varphi_i - \sigma^2/\mu - [...]$');
                end
                if any(strcmpi(OPTIONS, 'C'))
                    plot(SIGMA_LIST, phi_1_C_f(COMP_ID, :, p), 'b--', 'DisplayName', '$\varphi_i$');% = 2y_i\varphi_i - \sigma^2/\mu - [...]$');
                end

        end %NUM_POS
        myfigstyle(gcf, 12*num_col, 6*num_row, 10, 10);

        if any(strcmpi(OPTIONS, 'SAVE'))
            saveas(gcf,fullfile(SAVEPATH,['phi1','_comp',num2str(COMP_ID),'.fig']));
        end

        end %function
    end %methods
end


