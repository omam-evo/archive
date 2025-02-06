classdef Psi
    methods        
        %% Solve directly, no transformation
        function[psi_num] = get_psi_num(~, ras_tools, A, ALPHA, MU, LAM, TAU, sigma, Y_0)
            fprintf('> get_psi_num (sigma=%.3f) ...', sigma);
            VARTHETA = MU/LAM;
            %D_i = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_jni, sigma)));
            %E_x = @(x) ras_tools.expval_q_gain_fix_nonlin(A, ALPHA, x, Y_i, Y_jni, sigma);
            
            D_Q = @(s) sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_0, s))); 
            E_Q = @(s) sum(ras_tools.expval_q_gain(A, ALPHA, Y_0, s));
            
            p_ln = @(s)  1/sqrt(2*pi)/TAU/s * exp(-0.5*(log(s/sigma)/TAU)^2);
            %% Stat. properties
            m_transf = log(sigma);
            eval = exp(m_transf+0.5*TAU^2);
            mode = exp(m_transf-TAU^2);
            stddev = sqrt(exp(TAU^2-1)*exp(2*m_transf+TAU^2));
            
            s_min = max(mode-TAU,0);
            s_max = mode+TAU;

            PP_INV_q = @(s) E_Q(s) + D_Q(s)*norminv(VARTHETA);
            PP_q = @(s) normcdf((PP_INV_q(sigma) - E_Q(s))/D_Q(s));

            %% Debug
%             figure; hold on;
%             x_range = linspace(s_min, s_max, 1001);
%             p_ln_list = nan*x_range;
%             PP_q_list = nan*x_range;
%             terms_list = nan*x_range;
%             for i=1:length(x_range)
%                 x = x_range(i);
%                 p_ln_list(i) = p_ln(x);
%                 PP_q_list(i) = PP_q(x);
%                 terms_list(i) = (x-sigma)/sigma .* p_ln(x) .* PP_q(x);
%             end  
%             subplot(1,3,1); hold on; plot(x_range, p_ln_list);
%             subplot(1,3,2); hold on; plot(x_range, PP_q_list);
%             subplot(1,3,3); hold on; plot(x_range, terms_list);
%             dummy = 1;
            
            %% Solve
            psi_num = 1/VARTHETA * integral(@(s) (s-sigma)/sigma .* p_ln(s) .* PP_q(s), s_min, s_max, 'ArrayValued', 1);
            fprintf(' done.\n');      

        end   
        
        %% Integral transformed to be more robust w.r.t. small sigma
        % sigma-dep. moved into normcdf() -> maps to (0,1);
        % resolves integration issues for small sigma
        function[psi_num] = get_psi_num_v2(~, ras_tools, A, ALPHA, MU, LAM, TAU, sigma, Y_0)
            fprintf('> get_psi_num (sigma=%.3f) ...', sigma);
            VARTHETA = MU/LAM;
            %D_i = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_jni, sigma)));
            %E_x = @(x) ras_tools.expval_q_gain_fix_nonlin(A, ALPHA, x, Y_i, Y_jni, sigma);
            
            D_Q = @(s) sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_0, s))); 
            E_Q = @(s) sum(ras_tools.expval_q_gain(A, ALPHA, Y_0, s));
            
            p_ln = @(x)  1/sqrt(2*pi)/TAU * exp(-0.5*(x/TAU)^2);
            
            x_min = -6*TAU;
            x_max = 6*TAU;

            PP_INV_q = @(s) E_Q(s) + D_Q(s)*norminv(VARTHETA);
            PP_q = @(x) normcdf((PP_INV_q(sigma) - E_Q(sigma*exp(x)))/D_Q(sigma*exp(x)));

            %% Debug
%             figure; hold on;
%             x_range = linspace(x_min, x_max, 1001);
%             
%             p_ln_list = nan*x_range;
%             PP_q_list = nan*x_range;
%             terms_list = nan*x_range;
%             for i=1:length(x_range)
%                 x = x_range(i);
%                 p_ln_list(i) = p_ln(x);
%                 PP_q_list(i) = PP_q(x);
%                 terms_list(i) = (exp(x)-1) .* p_ln(x) .* PP_q(x);
%             end  
%             subplot(1,3,1); hold on; plot(x_range, p_ln_list);
%             subplot(1,3,2); hold on; plot(x_range, PP_q_list);
%             subplot(1,3,3); hold on; plot(x_range, terms_list);
%             dummy = 1;
            
            %% Solve
            psi_num = 1/VARTHETA * integral(@(x) (exp(x)-1) .* p_ln(x) .* PP_q(x), x_min, x_max, 'ArrayValued', 1);
            fprintf(' done.\n');     
        end   
        
        %% SMN C.30; SAR on sphere for small TAU in the limit N->inf
        function[psi_sph] = get_psi_sphere_signorm(~, MU, LAM, TAU, signorm)
            %FLAG_ZERO = 1; 
            %c_vartheta = e_vartheta_a_b(MU/LAM, 1, 0);
            %e_vt_11 = FLAG_ZERO*e_vartheta_a_b(MU/LAM, 1, 1);
            c_mu_lam = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
            e_mu_lam_11 = e_mu_lam_a_b_v2(MU, LAM, 1, 1);
            psi_sph = TAU^2*(0.5 + e_mu_lam_11 - c_mu_lam*signorm) ;
        end
        
        %% Analogous to SMN C.23; SAR on general function with E_Q, D_Q using large pop. identity giving e_vartheta_a_b instead of e_mu_lam_a_b
        function [psi,dDQ,DQ,dEQ] = get_psi_v1(~, ras_tools, A, ALPHA, TAU, sigma, Y, e_10, e_11)
      
            % --- Calculate derivative in each step ---
            %[D2Q, dD2Q] = ras_tools.variance_D2Q_dD2Q(A, ALPHA, Y, sigma);
            %[EQ, dEQ] = ras_tools.expect_EQ_dEQ(A, ALPHA, Y, sigma);
            
            % --- Pre-calculated derivative ---
            [~, dEQ] = ras_tools.EQ(A, ALPHA, sigma, Y);
            [D2Q, dD2Q] = ras_tools.D2Q(A, ALPHA, sigma, Y);
            
            DQ = sqrt(D2Q);
            dDQ = dD2Q/DQ/2;
            
            psi = TAU^2*(0.5 + e_11*sigma*dDQ/DQ - e_10*sigma*dEQ/DQ);  
            
        end
        
        %% Analogous to SMN C.23; SAR on general function with R-dependent EQ and D2Q
        function [psi,dDQ,DQ,dEQ] = get_psi_v1_R(~, ras_tools, A, ALPHA, TAU, sigma, R, N, e_10, e_11)
        
            % --- Pre-calculated derivative ---
            [~, dEQ] = ras_tools.EQ_Rdep(A, ALPHA, sigma, R, N);
            [D2Q, dD2Q] = ras_tools.D2Q_Rdep(A, ALPHA, sigma, R, N);
            
            DQ = sqrt(D2Q);
            dDQ = dD2Q/DQ/2;
            
            psi = TAU^2*(0.5 + e_11*sigma*dDQ/DQ - e_10*sigma*dEQ/DQ);  
            
        end
        
        %% V1 (first attempt)
        % bounds of s-integral not well argumented, only educated guess
        function[phi_1] = get_phi_via_lognormal_v1(~, ras_tools, A, ALPHA, MU, LAM, TAU, sigma, Y_0, Y_i, Y_jni)
            
            fprintf('> get_phi_via_lognormal (sigma=%.3f) ...', sigma);
            vartheta = MU/LAM;
            D_i = @(s) sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_jni, s)));
            E_i = @(s) sum(ras_tools.expval_q_gain(A, ALPHA, Y_jni, s));
            E_x = @(x,s) x^2 +2*x*Y_i - A*cos(ALPHA*Y_i)*cos(ALPHA*x) + A*sin(ALPHA*Y_i)*sin(ALPHA*x) + A*cos(ALPHA*Y_i) + E_i(s);
            
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_0, sigma))); 
            E_Q = sum(ras_tools.expval_q_gain(A, ALPHA, Y_0, sigma));
            
            p_x_s = @(x,s)  1/sqrt(2*pi)/s * exp(-(x/s).^2/2);
            p_ln = @(s)  1/sqrt(2*pi)/TAU/s * exp(-0.5*(log(s/sigma)/TAU)^2);
            
            % ---
            m_transf = log(sigma);
            %eval = exp(m_transf+0.5*TAU^2);
            mode = exp(m_transf-TAU^2);
            %stddev = sqrt(exp(TAU^2-1)*exp(2*m_transf+TAU^2));
            
            s_min = max(sigma - 6*exp(TAU),0); %max(mode-TAU,0);
            s_max = sigma + 6*exp(TAU);
%             figure;
%             fplot(p_ln, [s_min,s_max])
            
            PP_INV_q = E_Q + D_Q*norminv(vartheta);
            PP_q_x = @(x,s) normcdf((PP_INV_q - E_x(x,s))/D_i(s));
            
            phi_1 = -1/vartheta*integral(@(s) p_ln(s).*integral(@(x) x .* p_x_s(x,s) .* PP_q_x(x,s), -9*s, 9*s, 'ArrayValued', 1), s_min, s_max, 'ArrayValued', 1);
            fprintf(' done.\n');    
        end
     
        %% V1 with substitutions for rescaling      
        function[phi_1] = get_phi_via_lognormal_v2(~, ras_tools, A, ALPHA, MU, LAM, TAU, sigma, Y_0, Y_i, Y_jni)
            
            fprintf('> get_phi_via_lognormal (sigma=%.3f) ...', sigma);
            vartheta = MU/LAM;
            D_i = @(s) sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_jni, s)));
            E_i = @(s) sum(ras_tools.expval_q_gain(A, ALPHA, Y_jni, s));
            E_x = @(x,s) x^2 +2*x*Y_i - A*cos(ALPHA*Y_i)*cos(ALPHA*x) + A*sin(ALPHA*Y_i)*sin(ALPHA*x) + A*cos(ALPHA*Y_i) + E_i(s);
            
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_0, sigma))); 
            E_Q = sum(ras_tools.expval_q_gain(A, ALPHA, Y_0, sigma));
            
            p_z_y = @(z,y)  1/sqrt(2*pi) * exp(-0.5*(z/exp(y)).^2);
            p_ln = @(y)  1/sqrt(2*pi)/TAU * exp(-0.5*(y/TAU).^2);
            
            % ---
            y_min = -6*TAU;
            y_max = 6*TAU;
            
%             figure;
%             fplot(p_ln, [s_min,s_max])
            
            PP_INV_q = E_Q + D_Q*norminv(vartheta);
            PP_q_x = @(z,y) normcdf((PP_INV_q - E_x(sigma*z,sigma*exp(y)))/D_i(sigma*exp(y)));
            
            phi_1 = -sigma/vartheta * integral(@(y) exp(-y).*p_ln(y) .* integral(@(z) z .* p_z_y(z,y) .* PP_q_x(z,y), -6*exp(y), 6*exp(y), 'ArrayValued', 1), y_min, y_max, 'ArrayValued', 1);
            fprintf(' done.\n');    
        end
        
    end
end