% careful with definition of sqrt(..) in
% C_sqrt_inv = B'*diag(1./sqrt(diag(D)))*B;

function [PSA, p_t, delta_theta, norm2_p_dm, norm2_p_dS] = pc_psa_simpli(g, y_p, y_old, sigma, sigma_old, p_t, delta_theta, norm2_p_dm, norm2_p_dS, pc_input, cs, ds, E_chi, MU)
            N = size(y_p,1);

            dof = 2*N;
            beta = pc_input.beta_psa;
            mu_w = MU;
            c_m = 1; 
            alpha = pc_input.alpha_psa;
            
            %% Aki version with complex normalization
%             m1 = y_p;
%             m0 = y_old;
%             dm = m1-m0; 
%             dm = m1-m0;
%             dS = sigma^2*eye(N)-sigma_old^2*eye(N);  
%             gamma_t = 1;
%             gamma_t(g+1) = (1-beta)^2*gamma_t(g) + beta*(2-beta); 
%             E_I = N*c_m^2/mu_w + 2*N*(N-E_chi^2)/E_chi^2*gamma_t*(cs/ds)^2;
%             E_I = 1/mu_w; warning('Abgleich'); % Abgleich: ||s||^2=||p_m||^2
            
            %% zero if no CMA
            % E_I = E_I + 0.5*(1+8*gamma_t*(N-E_chi^2)/E_chi^2*(cs/ds)^2)* ...
            %    ((N^2+N)*cw^2/mu_w + (N^2+N)*cp*(2-cp)*c1*cw*mu_w*sum(w.^3)+ ...
            %    c1^2* (gamma_t^2*N^2 + N*(1-2*gamma_t+2*gamma_t^2) ) );

            %% PSA
            E_I = N/MU;
            dm_loc = (y_p-y_old) / sigma_old;
            dm_norm = dm_loc / sqrt(E_I);

            % dS_loc = 1/sqrt(2)*diag(SIG_SQRT_INV*dS*SIG_SQRT_INV);
            dS_loc = 1/sqrt(2)*(sigma^2/sigma_old^2-ones(N,1));
            dS_norm = dS_loc / sqrt(E_I);

            p_t(1:N,g+1) = (1-beta)*p_t(1:N,g) + sqrt(beta*(2-beta)) * dm_norm;
            p_t(N+1:dof,g+1) = (1-beta)*p_t(N+1:dof,g) + sqrt(beta*(2-beta)) * dS_norm;  
            delta_theta(:,g+1) = [dm_norm; dS_norm];           

            %% NA18, Eq.(15)
            norm2_p_dm(g+1) = norm(p_t(1:N,g+1)).^2;
            norm2_p_dS(g+1) = norm(p_t(N+1:dof,g+1)).^2;

            if norm2_p_dm(g+1)+norm2_p_dS(g+1) > alpha 
                PSA = 1;
            else
                PSA = -1;
            end
end