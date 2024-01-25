classdef IterVariants 
    methods
        %% Default full R-dep D_Q^2
        function [Var_q] = variance_full_R(~, A, ALPHA, sigma, R, N)
            sigma_Y = R/sqrt(N);
            Var_q = 2*N*sigma.^4 + 4*R.^2.*sigma.^2 + ...
                N*A^2/2 * (1-exp(-(ALPHA*sigma).^2)) .* (1-exp(-ALPHA^2*(sigma.^2 + 2*sigma_Y.^2))) + ...
                2*N*A*ALPHA^2*sigma.^2 .* exp(-0.5*ALPHA^2*(sigma.^2+sigma_Y.^2)) .* (sigma.^2 + 2*sigma_Y.^2);
        end
        %% Default full R-dep varphi_R
        function [phi_R, D_Q_R] = get_phi_R(obj, MU, A, ALPHA, sigma, R, N, c_vartheta)
            D_Q_R = sqrt(obj.variance_full_R(A, ALPHA, sigma, R, N));
            phi_R = c_vartheta * sigma.^2 .* (2*R.^2) ./ D_Q_R  .* (2 + ALPHA^2*A*exp(-0.5*ALPHA^2*(sigma.^2+R.^2/N))) - N*sigma.^2/MU;
        end
        
        %% Remove all exponentials
        function [phi_R, D_Q_R] = get_phi_R_v1(~, MU, A, ALPHA, sigma, R, N, c_vartheta)
            FLAG_ZERO = 0;
            sigma_Y = R/sqrt(N);
            Var_q = 2*N*sigma.^4 + 4*R.^2.*sigma.^2 + ...
                N*A^2/2 * (1-FLAG_ZERO*exp(-(ALPHA*sigma).^2)) .* (1-FLAG_ZERO*exp(-ALPHA^2*(sigma.^2 + 2*sigma_Y.^2))) + ...
                FLAG_ZERO*2*N*A*ALPHA^2*sigma.^2 .* exp(-0.5*ALPHA^2*(sigma.^2+sigma_Y.^2)) .* (sigma.^2 + 2*sigma_Y.^2);
            D_Q_R = sqrt(Var_q);
            
            phi_R = c_vartheta * sigma.^2 .* (2*R.^2) ./ D_Q_R  .* (2 + ...
                FLAG_ZERO*ALPHA^2*A*exp(-0.5*ALPHA^2*(sigma.^2+R.^2/N))) - N*sigma.^2/MU;
        end    
        
        %% Single exponential 
        function [phi_R, D_Q_R] = get_phi_R_v2(~, MU, A, ALPHA, sigma, R, N, c_vartheta)
            EXP = exp(-0.5*ALPHA^2*(sigma.^2+R.^2/N));
            
            Var_q = 2*N*sigma.^4 + 4*R.^2.*sigma.^2 + ...
                N*A^2/2 * (1-EXP);
            D_Q_R = sqrt(Var_q);
            
            phi_R = c_vartheta * sigma.^2 .* (2*R.^2) ./ D_Q_R  .* (2 + ALPHA^2*A*EXP) - N*sigma.^2/MU;
        end 
        
        %% Only transition (no sphere) 
        function [phi_R, D_Q_R] = get_phi_R_v3(~, MU, A, ALPHA, sigma, R, N, c_vartheta)
            EXP = exp(-0.5*ALPHA^2*(sigma.^2 + R.^2/N));
            ZERO = 0;
            
            sigma_Y = R/sqrt(N);
            Var_q = ZERO*(2*N*sigma.^4 + 4*R.^2.*sigma.^2) + ... 
                N*A^2/2 * (1-exp(-(ALPHA*sigma).^2)) .* (1-exp(-ALPHA^2*(sigma.^2 + 2*sigma_Y.^2)));
            D_Q_R = sqrt(Var_q);

            phi_R = c_vartheta * sigma.^2 .* (2*R.^2) ./ D_Q_R  .* (ZERO*2 + ALPHA^2*A*EXP) - N*sigma.^2/MU;
        end 
        
        %% Small R expansion for exp(..)
        function [phi_R, D_Q_R] = get_phi_R_v4(~, MU, A, ALPHA, sigma, R, N, c_vartheta)

            Var_q = (2+ALPHA^2*A)^2 * (sigma.^2 .* R.^2 + 0.5*N* sigma.^4 - ALPHA^4*A/(2+ALPHA^2*A)/2 * sigma.^4 .* R.^2);
            D_Q_R = sqrt(Var_q); % WARNING: will get complex for larger R
            
            EXP = 1 - 0.5*ALPHA^2 * (sigma.^2 + R.^2/N)  + 0.5*(0.5*ALPHA^2 * (sigma.^2 + R.^2/N)).^2;
            phi_R = c_vartheta * sigma.^2 .* (2*R.^2) ./ D_Q_R  .* (2 + ALPHA^2*A*EXP) - N*sigma.^2/MU;
        end
        
        %% D2_Q: neglect NA^2/2 term
        function [phi_R, D_Q_R] = get_phi_R_v5(~, MU, A, ALPHA, sigma, R, N, c_vartheta)
            sigma_Y = R/sqrt(N); 
            ZERO = 0;
            Var_q = 2*N*sigma.^4 + 4*R.^2.*sigma.^2 + ...
                ZERO*N*A^2/2 * (1-exp(-(ALPHA*sigma).^2)) .* (1-exp(-ALPHA^2*(sigma.^2 + 2*sigma_Y.^2))) + ...
                2*N*A*ALPHA^2*sigma.^2 .* exp(-0.5*ALPHA^2*(sigma.^2+sigma_Y.^2)) .* (sigma.^2 + 2*sigma_Y.^2);
            D_Q_R = sqrt(Var_q);
            
            phi_R = c_vartheta * sigma.^2 .* (2*R.^2) ./ D_Q_R  .* (2 + ALPHA^2*A*exp(-0.5*ALPHA^2*(sigma.^2+R.^2/N))) - N*sigma.^2/MU;
        end
        
        %% D2_Q: simplify NA^2/2 * [1-exp()]*[1-exp()]
        function [phi_R, D_Q_R] = get_phi_R_v6(~, MU, A, ALPHA, sigma, R, N, c_vartheta)
            sigma_Y = R/sqrt(N); 
            Var_q = 2*N*sigma.^4 + 4*R.^2.*sigma.^2 + ...
                N*A^2/2 * ( 1 - exp(-(ALPHA*sigma).^2) - exp(-ALPHA^2*(sigma.^2 + 2*sigma_Y.^2)) );
            D_Q_R = sqrt(Var_q);
            
            phi_R = c_vartheta * sigma.^2 .* (2*R.^2) ./ D_Q_R  .* (2 + ALPHA^2*A*exp(-0.5*ALPHA^2*(sigma.^2+R.^2/N))) - N*sigma.^2/MU;
        end   
        
        %% Local attractor using expected value
        function [phi_R_y, D_Q_R] = get_phi_local_attractor(obj, MU, A, ALPHA, sigma, R, N, c_vartheta, y_t, sigma_t)
            
            ev_sum = N*exp(-0.5*(ALPHA*sigma_t).^2).*(y_t*sin(ALPHA*y_t) + ALPHA*sigma_t.^2*cos(ALPHA*y_t));
            
            D_Q_R = sqrt(obj.variance_full_R(A, ALPHA, sigma, R, N));
            phi_R_y = c_vartheta * sigma.^2 ./ D_Q_R .*(4*R.^2 + ...
                exp(-0.5*ALPHA^2*sigma.^2)*2*ALPHA*A*ev_sum) + ...
                -N*sigma.^2/MU;
        end
        %% Local attractor using N*y*sin(alpha*y)
        function [phi_R_y, D_Q_R] = get_phi_local_attractor_v2(obj, MU, A, ALPHA, sigma, R, N, c_vartheta, y_t, sigma_t)
            
            sum_simplified = N*y_t*sin(ALPHA*y_t);
            
            D_Q_R = sqrt(obj.variance_full_R(A, ALPHA, sigma, R, N));
            phi_R_y = c_vartheta * sigma.^2 ./ D_Q_R .*(4*R.^2 + ...
                exp(-0.5*ALPHA^2*sigma.^2)*2*ALPHA*A*sum_simplified) + ...
                -N*sigma.^2/MU;
        end
        
    end
end