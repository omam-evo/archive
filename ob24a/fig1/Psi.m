classdef Psi
    methods 
        function [EV_q, d_EV_q] = expect_EQ_dEQ(~, A, ALPHA, Y, sigma)  
            
            syms s y
            f(s,y) = s^2 + A * cos(ALPHA*y)*(1 - exp(-(ALPHA*s)^2/2));           
            df = diff(f, s);
            f_fh = matlabFunction(f);   % Create function handle from symbolic expression => faster eval.
            df_fh = matlabFunction(df);
            EV_q = 0;
            d_EV_q = 0;
            for i=1:length(Y)
                EV_q = EV_q + f_fh(sigma,Y(i)); 
                d_EV_q = d_EV_q + df_fh(sigma,Y(i));
            end

        end
        function [Var_q, d_Var_q] = variance_D2Q_dD2Q(~, A, ALPHA, Y, sigma)  
            
            syms s y
            f(s,y) = 2*s^4 + 4*y^2*s^2 + ...
                A^2/2*(1 - exp(-(ALPHA*s)^2))*(1-cos(2*ALPHA*y)*exp(-(ALPHA*s)^2)) + ...
                2*A*ALPHA*s^2*exp(-0.5*(ALPHA*s)^2)*(ALPHA*s^2*cos(ALPHA*y) + 2*y*sin(ALPHA*y)); 
            df = diff(f, s);
            f_fh = matlabFunction(f);
            df_fh = matlabFunction(df);
            
            Var_q = 0;
            d_Var_q = 0;
            for i=1:length(Y)
                Var_q = Var_q + f_fh(sigma,Y(i)); 
                d_Var_q = d_Var_q + df_fh(sigma,Y(i));
            end

        end
        
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
        function[psi,dDQ,DQ,dEQ] = get_psi_v1(obj, A, ALPHA, MU, LAM, TAU, sigma, Y_0)
            
            fprintf('> get_psi_v1 (sigma=%.3f) ...', sigma);
            FLAG_ZERO = 0;
            
            %% Regular coefficients => compare with SMN result
            % c_vt = e_mu_lam_a_b_v2(MU, LAM, 1, 0);
            % e_vt_11 = e_mu_lam_a_b_v2(MU, LAM, 1, 1); 
            
            %% Asymtptotic coefficients
            c_vt = e_vartheta_a_b(MU/LAM, 1, 0);
            e_vt_11 = e_vartheta_a_b(MU/LAM, 1, 1);
            
            [D2Q, dD2Q] = obj.variance_D2Q_dD2Q(A, ALPHA, Y_0, sigma);
            [EQ, dEQ] = obj.expect_EQ_dEQ(A, ALPHA, Y_0, sigma);     
            DQ = sqrt(D2Q);
            dDQ = dD2Q/DQ/2;
            
            psi = TAU^2*(0.5 + e_vt_11*sigma*dDQ/DQ - c_vt*sigma*dEQ/DQ);
            fprintf(' done.\n');     
            
        end
    end
end