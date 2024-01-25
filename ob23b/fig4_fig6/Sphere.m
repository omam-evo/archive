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
                fplot(f, [0,1.2*x_2ndZero], 'k-')
                    plot(x_opt, f(x_opt), 'ro');
                    plot(x_2ndZero, f(x_2ndZero), 'go');
            end
        end % fct
       
        %% Solving: c*x/sqrt(1+x^2/(2*N)) - x^2/(2*mu) = 0
        % 
        function [x] = phiZero_medium_analytic(C_MU_LAM, MU, LAM, N, OPTIONS)
            if any(strcmpi(OPTIONS,'exact'))
                x = sqrt( sqrt(N*(8*C_MU_LAM^2*MU^2 + N)) - N);
            elseif any(strcmpi(OPTIONS,'approx1'))
                x = sqrt( sqrt(8)*C_MU_LAM*MU*sqrt(N) - N);
            end
        end
        
        %% Solving: d/dx c*x/sqrt(1+x^2/(2*N)) - x^2/(2*mu) = 0
        % 
        function [x] = signorm_opt_medium_analytic(C_MU_LAM, MU, LAM, N)
            sign1 = 1;
            sign2 = 1;
            S = ( sqrt(3)*(128*C_MU_LAM^6*MU^6*N^9 + 27*C_MU_LAM^4*MU^4*N^10)^(1/2) - 9*C_MU_LAM^2*MU^2*N^5 )^(1/3);
            %S = 9^(1/3)*(3*128/81)^(1/6)*C_MU_LAM*MU*N^(3/2); %% yields solution zero
            greenSqrt =  sqrt( 2*2^(1/3)/3^(2/3)*S - 8*2^(2/3)/3^(1/3)*C_MU_LAM^2*MU^2*N^3/S + N^2 );
            redSqrt = sqrt( -2*2^(1/3)/3^(2/3)*S + 8*2^(2/3)/3^(1/3)*C_MU_LAM^2*MU^2*N^3/S + ...
                2*N^3/greenSqrt + 2*N^2 );
            x = sqrt(sign1*0.5*greenSqrt +sign2*0.5*redSqrt -3*N/2 );
        end
           
    end % methods
end % class

