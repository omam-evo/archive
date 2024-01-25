%% Class containing the main results of Phi(1,2) Progress Rates for Rastrigin Function

classdef Phi_v2
    methods 
        
        %% Phi 1 via C_MU_LAM using D_Q for all i
        function[phi_1, ki, di] = get_phi_1_viaC_DQ(~, ras_tools, C_MU_LAM, A, ALPHA, sigma, Y, bool_di)
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y, sigma)));
            ki = 2*Y;
            if bool_di == 0
                %f_d1 = ras_tools.f_d1(Y, ALPHA, 0*A);
                di = 0;
            else
                di = ALPHA * A * sin(ALPHA*Y);
            end
            phi_1 = C_MU_LAM * sigma^2 / D_Q * (ki+di);
        end
        
        %% Phi 1 via C_MU_LAM using D_i for each i individually
        function[phi_1, k_i, d_i] = get_phi_1_viaC_Di(~, ras_tools, C_MU_LAM, A, ALPHA, sigma, Y, bool_di)
            [D2_k, f_p, k_i, d_i] = ras_tools.variance_q_gain_Df(A, ALPHA, Y, sigma, bool_di);
            D_k = sqrt(D2_k);
            phi_1 = C_MU_LAM * sigma^2 * f_p ./ D_k;

        end
        
        %% Phi 1 (N-dep) for sphere model; Eq. (6.54, TES)
        function[phi_1_R] = get_phi_1_sphere(~, MU, LAM, N, sigma_norm)
            C_MU_LAM = e_mu_lam_a_b_v2(MU,LAM,1,0);
            phi_1_R = C_MU_LAM * sigma_norm .* (1 + sigma_norm.^2/(2*MU*N))./sqrt(1+sigma_norm.^2/(2*N))./sqrt(1+sigma_norm.^2/(MU*N)) ...
                -N*(sqrt(1+sigma_norm.^2/(MU*N)) - 1);
        end
        
        %% Phi 1 according to A3
        function[phi_1, ki, exp_di] = get_phi_1(~, ras_tools, mu, lam, A, ALPHA, sigma, Y, bool_ki, bool_di)
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y, sigma)));
            vartheta = mu/lam;
            c_vartheta = e_vartheta_a_b(vartheta, 1, 0);
            ki = 2*Y;
            exp_di = exp(-(ALPHA*sigma)^2/2)*ALPHA*A*sin(ALPHA*Y);
            phi_1 = c_vartheta * sigma^2 / D_Q * (bool_ki*ki + bool_di*exp_di);
            % if any(imag(phi_1) ~= 0), pause, end
        end
        
        %% Phi 2 using Phi1(A3) and AR-approx. for E^(1,1) and E^(2)
        % => e_11, e_20 used to turn terms on/off
        function[phi_2] = get_phi_2(obj, ras_tools, mu, lam, A, ALPHA, sigma, Y, e_11, e_20)
            bool_ki = 1; 
            bool_di = 1;
            phi_1 = obj.get_phi_1(ras_tools, mu, lam, A, ALPHA, sigma, Y, bool_ki, bool_di);
            
            % --- Choose variance ---
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y, sigma)));
            % --- Variance D^2_k used before setting D_k=D_Q; D^2_k = ([ki+di]*sigma)^2+D_i^2 
            %bool_di = 0; 
            %[D2_k, ~, ~, ~] = ras_tools.variance_q_gain_Df(A, ALPHA, Y, sigma, bool_di);
            %D_Q = sqrt(D2_k);
            
            %% Get Phi2
            phi_2 = 2*Y.*phi_1 ...
                - sigma^2 / mu * (1 +(2*Y).^2 * sigma^2 ./ D_Q^2 * (e_11 + (mu-1)*e_20));
        end
        
        %% Phi 2 as function of R
        function[phi_2_R] = get_phi_2_R(obj, ras_tools, mu, A, alpha, sigma, R, N, c_vartheta)
            D_Q_R = sqrt(ras_tools.variance_q_gain_R(A, alpha, R, N, sigma));
            phi_2_R = c_vartheta*sigma^2/D_Q_R * 4*R^2* (1 + 0.5*alpha^2*A*exp(-0.5*alpha^2*(sigma^2+R^2/N))) - N*sigma^2/mu;
            
            % Test: calc phi with version 2 of eq.
            %sig_y = R/sqrt(N);
            %phi_v2 = c_vartheta*sigma^2/D_Q_R * (4*R^2 + 2*N*alpha^2*A*sig_y^2 * exp(-(alpha*sigma)^2/2)* exp(-(alpha*sig_y)^2/2) ) - ...
            %   N*sigma^2/mu;
        end
        
        %% TEST enviroment for phi_2 terms on/off etc.
        function[phi_2, T1_E2, T2_E2, T3_E11] = get_phi_2_lgPop(obj, ras_tools, mu, lam, A, ALPHA, sigma, Y_0, OPTIONS)
            
            vartheta = mu/lam;
            e_vt_10 = e_vartheta_a_b(vartheta, 1, 0); 
            e_vt_11 = e_vartheta_a_b(vartheta, 1, 1);
            e_vt_20 = e_vartheta_a_b(vartheta, 2, 0);
            
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_0, sigma)));
            expfac = exp(-0.5*(ALPHA*sigma)^2);
            
            %% V1 Terms from 1/mu^2*E^(2) and 1/mu^2*E^(1,1)
%             bool_exp = 1;
%             phi_1 = obj.get_phi_1(ras_tools, mu, lam, A, ALPHA, sigma, Y_0, bool_exp);
%             Term_x2_1 = sigma^2 / mu * (1 + e_vt_11*(2*Y_0).^2 * sigma^2 ./ D_Q^2);
%             Term_x2_2 = -sigma^2 / mu * e_vt_10/D_Q*(3*sigma^2 + A*cos(ALPHA*Y_0)*(1 - expfac + ALPHA^2*sigma^2*expfac) );
%             E2_overMuSq = Term_x2_1 + Term_x2_2;
%             
%             E11_overMuSq = (mu-1)/mu/2*phi_1.^2;       
%             
%             phi_2 = 2*Y_0.*phi_1 - E2_overMuSq - 2*E11_overMuSq;
           
            %% V2: calc. phi_1 and all terms as in report
            % TEST: 
            %   > setting T1=T2=T3=0 recovers sigma^2/mu
            %   > setting T1=T2=T3=0 matches B1/AR with sigma^2/mu + e_11 + e_20 up to e_vartheta
            E2_overMuSq = nan; E11_overMuSq = nan;
            k_i = 2*Y_0;
            d_i = ALPHA*A*sin(ALPHA*Y_0);
            phi_1 = e_vt_10 * sigma^2 / D_Q * (k_i + expfac*d_i);
            
            if any(strcmp(OPTIONS, 'L1'))
                ZERO = 0;
                T1_E2  = e_vt_11 * k_i.^2*sigma^2 / D_Q^2;
                T2_E2  = (-e_vt_10)/D_Q*(3*sigma^2 + A*cos(ALPHA*Y_0)*(1 - expfac + ALPHA^2*sigma^2*expfac));
                T3_E11  = (mu-1) * e_vt_20 * sigma^2 / D_Q^2 * (k_i + expfac*d_i).^2;
                phi_2 = 2*Y_0.*phi_1 - sigma^2/mu *(1 + T1_E2 + T2_E2 + T3_E11);
            elseif any(strcmp(OPTIONS, 'L2'))           
                phi_2 = 2*Y_0.*phi_1 - phi_1.^2;
            end
            
        end
        
        
        %% ========================= 
        % NUMERICAL SECTION
        
        function[phi_2_overNcomp] = phi_2_i_overNcomp(obj, ras_tools, A, ALPHA, MU, LAM, sigma, Y_0, OPTIONS)
            N = length(Y_0);
            phi_2_overNcomp = nan*Y_0;
            if any(strcmpi(OPTIONS, 'PARFOR'))
                parfor i=1:N
                    fprintf('> phi_2_i_overNcomp (i=%i) ...\n', i);
                    Y_i = Y_0(i); 
                    Y_jni = Y_0(i ~= 1:N);
                    [phi_2_overNcomp(i),~,~,~] = obj.get_phi_2_num_lgPop(ras_tools, A, ALPHA, MU, LAM, sigma, Y_0, Y_i, Y_jni);
                end
            else
                for i=1:N
                    fprintf('> phi_2_i_overNcomp (i=%i) ...\n', i);
                    Y_i = Y_0(i); 
                    Y_jni = Y_0(i ~= 1:N);
                    [phi_2_overNcomp(i),~,~,~] = obj.get_phi_2_num_lgPop(ras_tools, A, ALPHA, MU, LAM, sigma, Y_0, Y_i, Y_jni);
                end
            end
        end
        
        function[res] = get_phi_1_num(~, ras_tools, A, ALPHA, MU, LAM, sigma, Y_0, Y_i, Y_jni)
            
            fprintf('> get_phi_1_num (sigma=%.3f) ...\n', sigma);
            D_i = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_jni, sigma)));
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_0, sigma)));
            E_x = @(x) ras_tools.expval_q_gain_fix_nonlin(A, ALPHA, x, Y_i, Y_jni, sigma);
            E_Q = sum(ras_tools.expval_q_gain(A, ALPHA, Y_0, sigma));  
            
            % x INTEGRAL
            p_x = @(x)  1/sqrt(2*pi)/sigma * exp(-(x/sigma).^2/2);
            x_min = -6*sigma;
            x_max = 6*sigma;      
            PP_q_x = @(q, x) normcdf((q - E_x(x))/D_i);
            PP_INV_q = @(p) E_Q + D_Q*norminv(p);
            
            FAC = -(LAM-MU)*nchoosek(LAM,MU); %-LAM/beta(LAM-MU,MU)/MU
            t_min = 0.01;
            t_max = 0.99;
            
            intgrd_t = @(t) FAC * t.^(LAM-MU-1) .* (1-t).^(MU-1);
            intgrd_q_generic = @(t, x) PP_q_x(PP_INV_q(1-t), x);
            
            res = integral(@(x) x .* p_x(x) .* integral(@(t) intgrd_t(t).*intgrd_q_generic(t, x), t_min, t_max, 'ArrayValued', 1), x_min, x_max, 'ArrayValued', 1);
            fprintf(' done.\n');
        end
        
        function[phi_1_num_lgPop] = get_phi_1_num_lgPop(~, ras_tools, A, ALPHA, MU, LAM, sigma, Y_0, Y_i, Y_jni)
            fprintf('> get_phi_1_num_lgPop (sigma=%.3f) ...', sigma);
            vartheta = MU/LAM;
            D_i = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_jni, sigma)));
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_0, sigma)));
            E_x = @(x) ras_tools.expval_q_gain_fix_nonlin(A, ALPHA, x, Y_i, Y_jni, sigma);
            E_Q = sum(ras_tools.expval_q_gain(A, ALPHA, Y_0, sigma));
            
            p_x = @(x)  1/sqrt(2*pi)/sigma * exp(-(x/sigma).^2/2);
            x_min = -6*sigma;
            x_max = 6*sigma;
            
            PP_INV_q = @(p) E_Q + D_Q*norminv(p);
            PP_q_x = @(q, x) normcdf((q - E_x(x))/D_i);
            PINV_const = PP_INV_q(vartheta); % constant
            
            phi_1_num_lgPop = -1/vartheta*integral(@(x) x .* p_x(x) .* PP_q_x(PINV_const,x), x_min, x_max, 'ArrayValued', 1);
            fprintf(' done.\n');    
        end

        function[phi_2_num, phi_1_num_lgPop, E2_overMuSq, E11_overMuSq] = get_phi_2_num_lgPop(obj, ras_tools, A, ALPHA, MU, LAM, sigma, Y_0, Y_i, Y_jni)

            fprintf('> get_phi_2_num_lgPop (sigma=%.3f) ...\n', sigma);
            vartheta = MU/LAM;

            %% Single component
            D_i = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_jni, sigma)));
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_0, sigma)));
            E_x = @(x) ras_tools.expval_q_gain_fix_nonlin(A, ALPHA, x, Y_i, Y_jni, sigma);
            E_Q = sum(ras_tools.expval_q_gain(A, ALPHA, Y_0, sigma));

            %% Solution 1: -1/mu SUM(E[x_m;lam])= phi_1 
            %my_x =  @(x) x;
            %FAC = -1/vartheta/beta(lam-mu,mu);

            %% Solution 2: 1/mu^2 SUM(E[x^2_m;lam]) no approx. except normpdf() normcdf()
            my_x =  @(x) x.^2;
            FAC = 1/MU/vartheta/beta(LAM-MU,MU);       

            % x INTEGRAL
            p_x = @(x)  1/sqrt(2*pi)/sigma * exp(-(x/sigma).^2/2);
            x_min = -6*sigma;
            x_max = 6*sigma;      
            PP_q_x = @(q, x) normcdf((q - E_x(x))/D_i);
            PP_INV_q = @(p) E_Q + D_Q*norminv(p);

            % t INTEGRAL
            t_min = 0.001;
            t_max = 0.999;
            int_t = @(t) FAC * t.^(LAM-MU-1) .* (1-t).^(MU-1);
            int_ft = @(t, x) PP_q_x(PP_INV_q(1-t), x);

            %% Solve: 1/mu^2 SUM(E[x^2_m;lam])
            fprintf('\t > numeric E2_overMuSq ...', sigma);
                E2_overMuSq = ... 
                    integral(@(x) my_x(x) .* p_x(x) .* ... 
                        integral(@(t) int_t(t).*int_ft(t, x), t_min, t_max, 'ArrayValued', 1), ...
                    x_min, x_max, 'ArrayValued', 1);
            fprintf(' done.\n');

            %% Solve: \varphi_i by large Pop. Approx.
%             PINV_const = PP_INV_q(vartheta); % constant
%             fprintf('\t numeric phi_i [lgPop] ...', sigma);
%             phi_1_num_lgPop = -1/vartheta*integral(@(x) x .* p_x(x) .* PP_q_x(PINV_const,x), x_min, x_max, 'ArrayValued', 1);
%             fprintf(' done.\n');    
            fprintf('\t ');
            phi_1_num_lgPop = obj.get_phi_1_num_lgPop(ras_tools, A, ALPHA, MU, LAM, sigma, Y_0, Y_i, Y_jni);

            E11_overMuSq = (MU-1)/MU/2*phi_1_num_lgPop^2;

            phi_2_num = 2*Y_i*phi_1_num_lgPop - E2_overMuSq - 2*E11_overMuSq;
        end

        %% EXPERIMENTAL
        
        function[phi_2_mean, phi_2_stderr ] = phi_2_experiment(obj, TRIALS, FIT, MU, LAM, Y_0, SIGMA_0)

        N = length(Y_0);
        N_G = 1; SIGMA_STOP = nan; F_STOP = nan; TAU = 0; Y_HAT=0*Y_0;
        [phi_1_trial, phi_2_trial] = deal(zeros(N, TRIALS), zeros(N, TRIALS));

        parfor t = 1:TRIALS

            rstream = RandStream('threefry4x64_20');
            rstream.Substream = t;

            %F_0 = FIT(Y_0);
            %R_0 = norm(Y_0 - Y_HAT);

            [y_opt, F_opt_vec, delta_r_vec, sigma_vec, gen] = muComLam_sigSA_randStream(rstream, FIT, N, MU, LAM, Y_0, Y_HAT, SIGMA_0, TAU, SIGMA_STOP, F_STOP, N_G, 0);

            % Phi_1, Phi_2, Phi_R, Q-Gain
            % phi_1_trial(:, t) = Y_0(:) - y_opt(:);
            phi_2_trial(:, t) = Y_0(:).^2 - y_opt(:).^2;
            %phi_R2_trial(1, t) = R_0^2 - norm(y_opt - Y_HAT)^2;
            
        end % trials
        
        % Phi 1
        %phi_1_mean(:, s, p) = mean(phi_1_trial, 2);
        %phi_1_stderr(:, s, p) = std(phi_1_trial, 0, 2)/sqrt(TRIALS); 
        % Phi 2
        phi_2_mean = mean(phi_2_trial, 2);
        phi_2_stderr = std(phi_2_trial, 0, 2)/sqrt(TRIALS);  
        % Phi R
        %phi_R2_mean(s, p) = mean(phi_R2_trial);
        %phi_R2_stderr(s, p) = std(phi_R2_trial, 0, 2)/sqrt(TRIALS);  
 
    end 
        
    end
end