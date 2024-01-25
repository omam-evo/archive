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
        function[phi_1, ki, di] = get_phi_1_viaC_Di(~, ras_tools, C_MU_LAM, A, ALPHA, sigma, Y, bool_di)
            N = length(Y);
            phi_1 = nan*Y;
            
            for i=1:N
                Y_i = Y(i);
                Y_jni = Y(i ~= 1:N);
                D_i = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y_jni, sigma)));
                ki = 2*Y_i;
                if bool_di == 0
                    di = 0;
                else
                    di = ALPHA * A * sin(ALPHA*Y_i);
                end
                f_p = ki+di;
                D_k = sqrt((f_p*sigma)^2 + D_i^2);
                phi_1(i) = C_MU_LAM * sigma^2 / D_k * f_p;
            end

        end
        %% Phi 1 according to A3
        function[phi_1, ki, exp_di] = get_phi_1(~, ras_tools, mu, lam, A, ALPHA, sigma, Y, bool_exp)
            D_Q = sqrt(sum(ras_tools.variance_q_gain_simplified(A, ALPHA, Y, sigma)));
            theta = mu/lam;
            c_theta = 1/sqrt(2*pi)/theta * exp(-1/2*(norminv(theta))^2);
            ki = 2*Y;
            exp_di = bool_exp*exp(-(ALPHA*sigma)^2/2)*ALPHA*A*sin(ALPHA*Y);
            phi_1 = c_theta * sigma^2 / D_Q * (ki + exp_di);
%             if any(imag(phi_1) ~= 0)
%                 pause
%             end
        end
        function[phi_2] = get_phi_2(obj, ras_tools, mu, lam, A, ALPHA, sigma, Y, bool_exp)
            phi_1 = obj.get_phi_1(ras_tools, mu, lam, A, ALPHA, sigma, Y, bool_exp);
            phi_2 = 2*Y.*phi_1 - sigma^2 / mu;
        end
    end
end