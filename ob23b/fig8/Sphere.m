classdef Sphere 
    methods     (Static=true)
        function [res] = phi(sigma_norm, C_MU_LAM, MU, LAM, N, OPTIONS)

            syms x
            if any(strcmpi(OPTIONS,'full'))
                f(x) = C_MU_LAM * x * (1 + x^2/(2*MU*N))/sqrt(1+x^2/(2*N))/sqrt(1+x^2/(MU*N)) - N*(sqrt(1+x^2/(MU*N)) - 1);
            elseif any(strcmpi(OPTIONS,'medium'))
                f(x) = C_MU_LAM * x / sqrt(1+x^2/(2*N)) - x^2/(2*MU);
            else
                warning('Eq. type not defined.')
            end
            
            res = double(f(sigma_norm));
            
        end % fct
        
        function [x_opt, x_2ndZero] = signorm_opt_sph(C_MU_LAM, MU, LAM, N, OPTIONS)

            %C_MU_LAM = e_mu_lam_a_b_v2(MU,LAM,1,0);
            syms x
            if any(strcmpi(OPTIONS,'full'))
                f(x) = C_MU_LAM * x * (1 + x^2/(2*MU*N))/sqrt(1+x^2/(2*N))/sqrt(1+x^2/(MU*N)) - N*(sqrt(1+x^2/(MU*N)) - 1);
            elseif any(strcmpi(OPTIONS,'medium'))
                f(x) = C_MU_LAM * x / sqrt(1+x^2/(2*N)) - x^2/(2*MU);
            else
                warning('Eq. type not defined.')
            end
            df = diff(f, x);
            assume(x>0);
            x_opt = double(solve(df == 0, x));

            assume(x>x_opt);
            x_2ndZero = double(solve(f == 0, x));

            if any(strcmpi(OPTIONS,'PLOT'))
                figure; hold on;
                fplot(f, [0,20], 'k-')
                    plot(x_opt, f(x_opt), 'ro');
                    plot(x_2ndZero, f(x_2ndZero), 'go');
            end
        end % fct
        
    end % methods
end % class

