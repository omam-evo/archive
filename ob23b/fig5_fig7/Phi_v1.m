%% Working on Phi derivations incl. numeric solutions and different analytic approximation

classdef Phi
    properties
        MU; LAM; SIGMA; MODE_E_x;
        fp_i; fp_all; k_i; d_i; D_f; D_k;
        E_i; D_i; E_Q; D_Q; D_Q_smSig; D_Q_lgSig;
        m; s;
        x_lim; E_x; I_x;
        PP_Q; PP_q_x; PP_INV_q; 
        phi_genPQ_num;  % No approx. besides P(q)=normcdf(q-Mq/Sq) => see documentation when exact/approx
        phi_lgPop_num;  % Large Pop.Approx with hyperg. Fct.
        phi_lgPop;      % Large Pop.Approx; Linearized i-th component of QGain
        phi_c;          % Linearized i-th component of QGain; D_i~D_Q for large N; q-integration over (-inf, inf)
        phi_c_num;      % No approx. besides P(q)=normcdf(q-Mq/Sq) => see documentation when exact/approx
        phi_c_smSig;
        phi_c_lgSig;
    end
    
    methods   
        function [obj] = Phi(MU, LAM, ras_tools, A, ALPHA, Y, Y_i, Y_inj, sigma, MODE_E_x)
            
            obj.MU = MU;
            obj.LAM = LAM;
            obj.SIGMA = sigma;
            obj.MODE_E_x = MODE_E_x;
            
            obj.x_lim = sigma*6;
            obj.I_x = @(x)  1/sqrt(2*pi*sigma^2) * x .* exp(-(x/sigma).^2/2);
            
            if strcmpi(MODE_E_x, 'linear')
                obj.E_x = @(x) ras_tools.expval_q_gain_fix(A, ALPHA, x, Y_i, Y_inj, sigma, 1);
            elseif strcmpi(MODE_E_x, 'nonlinear')
                obj.E_x = @(x) ras_tools.expval_q_gain_fix_nonlin(A, ALPHA, x, Y_i, Y_inj, sigma);
            elseif strcmpi(MODE_E_x, 'quadr')
                obj.E_x = @(x) ras_tools.expval_q_gain_fix(A, ALPHA, x, Y_i, Y_inj, sigma, 2);
            elseif strcmpi(MODE_E_x, 'testing')
                obj.E_x = @(x) ras_tools.expval_q_gain_fix_testing(A, ALPHA, x, Y_i, Y_inj, sigma);
            end
            
            obj.E_i = sum(ras_tools.expval_q_gain(A, ALPHA, Y_inj, sigma));
            obj.E_Q = sum(ras_tools.expval_q_gain(A, ALPHA, Y, sigma));

            obj.D_i = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_inj, sigma)));
            obj.D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y, sigma)));

            obj.PP_Q = @(q) normcdf((q - obj.E_Q) / obj.D_Q);
            obj.PP_q_x = @(q, x) normcdf((q - obj.E_x(x))/obj.D_i);
            obj.PP_INV_q = @(p) obj.E_Q + obj.D_Q*norminv(p);
            
            obj.fp_i = ras_tools.f_d1(Y_i, ALPHA, A);
            obj.fp_all = ras_tools.f_d1(Y, ALPHA, A);
            obj.k_i = 2*Y_i;
            obj.d_i = ALPHA*A*sin(ALPHA*Y_i);
            obj.D_f = sqrt(obj.D_i^2 + (obj.fp_i*sigma)^2);
            obj.D_k = sqrt(obj.D_i^2 + (obj.k_i*sigma)^2);
            
            %% Approx.
            obj.D_Q_smSig = sqrt(sum(ras_tools.variance_q_gain_smallSig(A, ALPHA, Y, sigma, 0) ));
            obj.D_Q_lgSig = sqrt(sum(ras_tools.variance_q_gain_largeSig(A, ALPHA, Y, sigma)));
            obj.m = (obj.E_Q-obj.E_i + obj.D_Q*norminv(obj.MU/obj.LAM)) * (obj.k_i*obj.SIGMA^2) / obj.D_k^2;
            obj.s = obj.SIGMA*obj.D_i/obj.D_k;
            
        end
        
        %% Numeric solution of phi_1 via progress coeff. (I1 in Tech. Report)
        % MODE_E_x = 'nonlinear' and SET_E_ALL_APPROX == 0
        %    => EXACT numeric solution as get_phi_genericPQ_num()
        % MODE_E_x = 'linear' and SET_E_ALL_APPROX == 1
        %    => Numeric solution & validation of get_phi_viaC() => no deviations detected
        function[res, obj] = get_phi_SumOverM_num(obj, SET_E_ALL_APPROX)         
            %% Mutation integrand         
            if SET_E_ALL_APPROX == 0
                % Integration performed  
                x_limit = 6 * obj.SIGMA;
                p_q_x = @(x, q) 1/sqrt(2*pi*obj.D_i^2) * exp(-(q - obj.E_x(x)).^2 / (obj.D_i^2) / 2);
                intgrd_x = @(x, q) obj.I_x(x).*p_q_x(x, q);
                I_q_1 = @(q) integral(@(x) intgrd_x(x, q), -x_limit, x_limit, 'ArrayValued', 1);
            else
                % Use already calculated integral to validate analytical derivation => Insert E_i = E_Q and D_i = D_Q
                % Tech.Rep. ~(93) before "The integration range of q was set to..."
                I_q_1 = @(q) 1/sqrt(2*pi)* obj.fp_i * obj.SIGMA^2 / obj.D_Q^2 * ((q - obj.E_Q)/obj.D_Q) .* exp(-1/2*((q - obj.E_Q)/obj.D_Q).^2);
            end
            
            %% Population dependent integrand
            I_q_2 = @(q, m) (-obj.LAM/obj.MU)*factorial(obj.LAM-1)/factorial(m-1)/factorial(obj.LAM-m) * obj.PP_Q(q).^(m-1) .* (1-obj.PP_Q(q)).^(obj.LAM-m);
            
            %% Combining integrands
            intgrd_q = @(q, m) I_q_1(q).*I_q_2(q, m);
            q_lim = obj.D_Q*10;
            q_lim_min = -q_lim + obj.E_Q;
            q_lim_max =  q_lim + obj.E_Q;
            
            fprintf('\t get_phi_SumOverM_num  (sigma=%.3f) ...', obj.SIGMA);
            res = 0;
            for m=1:obj.MU
                res = res + integral(@(q) intgrd_q(q, m), q_lim_min, q_lim_max);
            end
            obj.phi_c_num = res;
            fprintf(' done.\n');
        end
        
        %% Numeric solution of phi_1 via generic P_Q(q) integration (I3 in Tech. Report)
        % MODE_E_x = 'nonlinear': exact solution as get_phi_SumOverM_num()
        % NO (!) large pop. approx. applied
        function[res, obj] = get_phi_genericPQ_num(obj)
            
            FAC = -(obj.LAM-obj.MU)*nchoosek(obj.LAM,obj.MU);
            t_lim_min = 0;
            t_lim_max = 1;
            
            intgrd_t = @(t) FAC * t.^(obj.LAM-obj.MU-1) .* (1-t).^(obj.MU-1);
            intgrd_q_generic = @(t, x) obj.PP_q_x(obj.PP_INV_q(1-t), x);
            
            fprintf('\t get_phi_genericPQ_num (sigma=%.3f) ...', obj.SIGMA);
                res = integral(@(x) obj.I_x(x) .* integral(@(t) intgrd_t(t).*intgrd_q_generic(t, x), t_lim_min, t_lim_max, 'ArrayValued', 1), -obj.x_lim, obj.x_lim, 'ArrayValued', 1);
                obj.phi_genPQ_num = res;
            fprintf(' done.\n');
        end
        
        %% Numeric solution of phi_1 with large Pop. approx. (I4 in Tech. Report)
        % MODE_E_x = 'linear': numeric solution & validation of get_phi_large_pop() => no deviations detected
        function[res, obj] = get_phi_large_pop_num(obj)
            THETA = obj.MU/obj.LAM;
            v_hat = 1 - THETA; %asymt. correct
            intgrd_q_generic_vHat = @(x) obj.PP_q_x(obj.PP_INV_q(1-v_hat), x);
            fprintf('\t get_phi_large_pop_num (sigma=%.3f) ...', obj.SIGMA);
                res = -(1/THETA) * integral(@(x) obj.I_x(x) .* intgrd_q_generic_vHat(x), -obj.x_lim, obj.x_lim, 'ArrayValued', 1);
                obj.phi_lgPop_num = res;
            fprintf(' done.\n');
        end
       
       %% phi_1 via large pop. approx. (linear Taylor term for Q)
       function[res, obj] = get_phi_large_pop(obj)  
            E_Qi = 1*(obj.E_Q-obj.E_i);   %turn on/off i-th component (minor influence)
            sig_tilde = obj.D_f; %obj.D_f,obj.D_Q; SETTING: sig_tilde = obj.D_Q, E_Qi = 0*(obj.E_Q-obj.E_i) => phi_1 via C
            THETA = obj.MU/obj.LAM;
            res = 1/sqrt(2*pi)/THETA * (obj.fp_i*obj.SIGMA^2) / sig_tilde * exp(-1/2*(E_Qi+obj.D_Q*norminv(THETA))^2/sig_tilde^2);
            obj.phi_lgPop = res;
       end
       
       %% phi_1 via progress coeff. (linear Taylor term for Q)
       function[res_c, obj] = get_phi_viaC(obj, C_MU_LAM, mode)
            if strcmpi(mode,'fpi')
                res_c = C_MU_LAM * obj.SIGMA^2 / obj.D_Q * obj.fp_i;
            elseif strcmpi(mode,'ki')
                res_c = C_MU_LAM * obj.SIGMA^2 / obj.D_Q * obj.k_i;
            else %undefined
                assert(false);
            end
            obj.phi_c = res_c;
            
            %res_smSig = C_MU_LAM * obj.SIGMA * obj.fp_i/norm(obj.fp_all);
            %obj.phi_c_smSig = res_smSig;        
            %res_lgSig = C_MU_LAM * (obj.k_i*obj.SIGMA^2)/(obj.D_Q + 0*obj.D_Q_lgSig);
            %obj.phi_c_lgSig = res_lgSig;
       end
       
        %% Numeric solution of phi_1: Phi(A+B)
        % Linear part yields same as get_phi_large_pop(), if obj.k_i set to obj.fp_i (not 2*Y_i)
        function[res, obj] = get_phi_PhiAB_num(obj, A, ALPHA, Y_i)
            THETA = obj.MU/obj.LAM;
            
            %% Linear
            E_Qi = 1*(obj.E_Q-obj.E_i); 
            c1 = -obj.k_i/obj.D_i;
            c0 = (E_Qi + obj.D_Q*norminv(THETA)) /obj.D_i;
            g = @(x) c1 * x + c0;
            Phi_lin = @(x) normcdf(g(x));
            fprintf('\t get_phi_PhiExpand_num (0th) (sigma=%.3f) ...', obj.SIGMA);
                res_lin = -(1/THETA) * integral(@(x) obj.I_x(x) .* Phi_lin(x), -obj.x_lim, obj.x_lim, 'ArrayValued', 1);
                res = res_lin;
            fprintf(' done.\n');
            
            %% Non-Linear
            a_1 = 1;
            c_i = 1*A*cos(ALPHA*Y_i);
            s_i = 1*A*sin(ALPHA*Y_i);
            h = @(x) (-1)*(a_1*x.^2 + s_i*sin(ALPHA*x) + c_i*(1-cos(ALPHA*x)))/obj.D_i;
            Integrand_NonLin = @(x) h(x) .* exp(-1/2*g(x).^2)/sqrt(2*pi);
            
            % s = obj.SIGMA*obj.D_i/obj.D_k;
            % m = (E_Qi + obj.D_Q*norminv(THETA)) /obj.D_k * (obj.k_i*obj.SIGMA^2)/obj.D_k;
            x_min = obj.m - 6*obj.s;
            x_max = obj.m + 6*obj.s;
            
            fprintf('\t get_phi_PhiExpand_num (1st) (sigma=%.3f) ...', obj.SIGMA);
                res_nonlin = -(1/THETA) * integral(@(x) obj.I_x(x) .* Integrand_NonLin(x), x_min, x_max, 'ArrayValued', 1);
                res = res + res_nonlin;
            fprintf(' done.\n');
            
%             Integrand_NonLin_2nd = @(x) -1/(2*sqrt(2*pi)) * g(x).*(h(x).^2).*exp(-1/2*g(x).^2);
%             fprintf('\t get_phi_PhiExpand_num (2nd) (sigma=%.3f) ...', obj.SIGMA);
%                 res_nonlin_2nd = -(1/THETA) * integral(@(x) obj.I_x(x) .* Integrand_NonLin_2nd(x), x_min, x_max, 'ArrayValued', 1);
%                 res = res + res_nonlin_2nd;
%             fprintf(' done.\n');
        end
        function[res, obj] = get_phi_PhiAB(obj, ras_tools, A, ALPHA, Y_i)
            %% General definitions
            THETA = obj.MU/obj.LAM;
            E_Qi = 1*(obj.E_Q-obj.E_i);  
            a_1 = 1;
            c_i = 1*A*cos(ALPHA*Y_i);
            s_i = 1*A*sin(ALPHA*Y_i);          
            
            %% Linear
            I_A = 1/sqrt(2*pi)/THETA * (obj.k_i*obj.SIGMA^2) / obj.D_k * exp(-1/2*(E_Qi+obj.D_Q*norminv(THETA))^2/obj.D_k^2);
            
            %% Version 1: using the results from quadratic completion
%                 c0 = (E_Qi + obj.D_Q*norminv(THETA)) / obj.D_i;
%                 c1 = -obj.k_i / obj.D_i;
%                 beta =  1/obj.SIGMA^2 + c1^2;
%                 m = -c0*c1/beta;
%                 s = 1/sqrt(beta);
%                 C = exp((c0*c1)^2/2/beta - c0^2/2);
%                 T1 = C*s/sqrt(2*pi)/THETA/obj.SIGMA/obj.D_i;
%                 T2 = a_1*ras_tools.EV_x3(m, s) + s_i*ras_tools.EV_x1sin1(ALPHA, m, s) + c_i*(m - ras_tools.EV_x1cos1(ALPHA, m, s));
                
            %% Version 2: exact intermediate result 
                s_temp = obj.SIGMA*obj.D_i/obj.D_k;
                m_temp = (E_Qi + obj.D_Q*norminv(THETA)) * (obj.k_i*obj.SIGMA^2) / obj.D_k^2;
                C = exp(1/2*(E_Qi + obj.D_Q*norminv(THETA))^2/obj.D_i^2 * ((obj.k_i*obj.SIGMA/obj.D_k)^2-1) );
                T1 = 1/sqrt(2*pi)/THETA/obj.D_k * C;
                T2 = a_1*ras_tools.EV_x3(m_temp, s_temp) + s_i*ras_tools.EV_x1sin1(ALPHA, m_temp, s_temp) + c_i*(m_temp - ras_tools.EV_x1cos1(ALPHA, m_temp, s_temp));

            %% Version 3: for EValues O(sigma^3) neglected; validated using numerical calc. and a_1=0, c_i=0, E_Qi=0
%                 T1 = 1/sqrt(2*pi)/THETA/obj.D_k * C;
%                 T2 = s_i*exp(-1/2*(ALPHA*s)^2)*(m*sin(ALPHA*m)+ALPHA*s^2*cos(ALPHA*m)); %s_i*ras_tools.EV_x1sin1(ALPHA, m, s)
            
            %% Result
            I_B = T1 * T2;
            
            %% Res
            res = I_A + I_B;
        end
        
        function[res, I_A, I_B, obj] = get_phi_PhiAB_approx(obj, ras_tools, A, ALPHA, Y_i)
            %% General definitions
            THETA = obj.MU/obj.LAM;
            a_1 = 1;
            c_i = 1*A*cos(ALPHA*Y_i);
            s_i = 1*A*sin(ALPHA*Y_i); 
            m_bool = 0;
            
            %%
            % Terms simplified by setting a_1=0, c_i=0, E_Qi=0, D_Q=D_i=D_k
            % One of the variances D_* must be chosen:
            E_Qi = 0*(obj.E_Q-obj.E_i);  
            D_Q_temp = obj.D_Q_lgSig; %obj.D_Q; obj.D_Q_smSig
            D_i_temp = obj.D_Q_lgSig; %obj.D_i;
            D_k_quad_temp = obj.D_Q_lgSig; %obj.D_k;
            
            s_temp = obj.SIGMA * D_i_temp/D_k_quad_temp;
            m_temp = m_bool*(E_Qi + D_Q_temp*norminv(THETA)) * (obj.k_i*obj.SIGMA^2) / D_k_quad_temp^2;
      
            % common coefficient
            c_theta = 1/sqrt(2*pi)/THETA * exp(-1/2*(E_Qi+D_Q_temp*norminv(THETA))^2/D_k_quad_temp^2);
            
            %% Linear
            I_A = c_theta * obj.k_i * obj.SIGMA^2 / D_k_quad_temp;
            I_B = c_theta * (1/D_k_quad_temp) * (a_1*ras_tools.EV_x3(m_temp, s_temp) + s_i*ras_tools.EV_x1sin1(ALPHA, m_temp, s_temp) + c_i*(m_temp - ras_tools.EV_x1cos1(ALPHA, m_temp, s_temp)));
            res = I_A + I_B;
           
        end
        function[res, T1, T2, obj] = get_phi_PhiAB_A3(obj, A, ALPHA, Y_i)
            THETA = obj.MU/obj.LAM;
            c_theta = 1/sqrt(2*pi)/THETA * exp(-1/2*(norminv(THETA))^2);
            T1 = c_theta * obj.SIGMA^2 / obj.D_Q * 2*Y_i;
            T2 = c_theta * obj.SIGMA^2 / obj.D_Q * exp(-(ALPHA*obj.SIGMA)^2/2)*ALPHA*A*sin(ALPHA*Y_i);
            res = T1 + T2;
%             res = c_theta * obj.SIGMA^2 / obj.D_Q * (2*Y_i + exp(-(ALPHA*obj.SIGMA)^2/2)*ALPHA*A*sin(ALPHA*Y_i));
%             warning('get_phi_PhiAB_A3(): D_k_quad')
        end
        
        function[res] = get_phi2_B1(obj, phi_1, Y_i, e_11, e_20)
            T_1 = 2*Y_i*phi_1;
            T_2 = -obj.SIGMA^2 / obj.MU;
            T_3 = -(2*Y_i)^2*obj.SIGMA^4/obj.MU/obj.D_k^2 * (e_11 + (obj.MU-1)*e_20);
            res = T_1 + T_2 + T_3;
        end 
        
        function[res] = get_phi2(obj, phi_1, Y_i)
            T_gain = 2*Y_i*phi_1;
            T_loss = -obj.SIGMA^2 / obj.MU;
            res = T_gain + T_loss;
        end
        
    end % methods
end % class