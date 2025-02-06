%% CSA to call externally
% C = 1/sqrt(N_DIM);
% D = sqrt(N_DIM);
% E_chi = sqrt(N_DIM)*(1 - 1/(4*N_DIM) + 1/(21*N_DIM^2));

function [y_rec, f_g, r_g, sigma_g, y_g, pc_data] = csa_es_meas(pc_input, FIT, CSA_MODE, N, MU, LAM, Y_0, SIGMA_0, S_0, C_handle, D_handle, E_chi, SIGMA_STOP, R_STOP, F_STOP, FEVAL_STOP, G, verbose)

    f_g = nan*zeros(G,1);
    sigma_g = nan*ones(G,1);
    r_g = nan*ones(G,1);
    %y_g = nan*ones(N,G); 
    y_g = nan; % assign dummy value as output
    s_g = nan*ones(N,G); 
    z_g = nan*ones(N,G); 
    mu_g = nan*ones(G,1);
    THETA = MU/LAM;
    feval_budget = FEVAL_STOP;

    y_rec = Y_0;
    s = S_0*ones(N,1);
    sigma = SIGMA_0;
    
    % --- needed for non-vectorized version ---
    %z_tilde = nan*ones(N, LAM); 
    %y_tilde = nan*ones(N, LAM); 
    %f_tilde = nan*ones(1, LAM);

    %% Measures 
    D2 = nan*zeros(G,1); %D2(1)=0; % INFNOISE
    Q = nan*zeros(G,1); Q(1) = 0; % QGAIN
    dF = nan*zeros(G,1); dF(1) = 0; % ROSAES
    Z = nan*zeros(N,G); Z(:,1) = 0; % LENGTHZ
    P_LENGTHZ = nan*zeros(G, 1); % LENGTHZ
    N_up = nan*zeros(G, 1); % APOP
    H = nan*zeros(G, 1); % HEMI
    PH = nan*zeros(G, 1); % HEMI
    P_INFNOISE = nan*zeros(G, 1); % INFNOISE
    P_RANK = nan*zeros(G, 1); % RANK
    
    %% F-values & gen
    F_mean_g = nan*zeros(G,1); 
    F_rec_g = nan*zeros(G, 1); 
    F_med_g = nan*zeros(G, 1);
    wait = pc_input.wait;

    %% NA18
    dof = 2*N;
    p_t = nan*zeros(dof, G);
    p_t(:,1) = 0;
    norm2_p_dm = nan*zeros(G,1);  
    norm2_p_dS = nan*zeros(G,1); 
    delta_theta = nan*zeros(dof, G);  
    
    %% Init
    f_g(1) = FIT.get_f(y_rec);
    [F_rec_g(1), F_med_g(1),  F_mean_g(1)] = deal(f_g(1),f_g(1),f_g(1));
    sigma_g(1) = sigma;
    r_g(1) = norm(y_rec-FIT.y_hat);
    %y_g(:,1) = y;
    mu_g(1) = MU;
    s_g(:,1) = s;

    for g = 1:G-1

        z_tilde = randn(N, LAM);
        y_tilde = repmat(y_rec, 1, LAM) + sigma*z_tilde;
        f_tilde = FIT.get_f(y_tilde, g, norm(y_rec-FIT.y_hat));
        feval_budget = feval_budget-LAM;
        
        % for i=1:LAM
        %     z_tilde(:, i) = randn(N, 1);
        %     y_tilde(:, i) = y_p + sigma*z_tilde(:, i);
        %     f_tilde(i) = FIT(y_tilde(:,i));
        % end

        [~, idx] =  sort(f_tilde, 'ascend');
        y_old = y_rec;
        y_rec = mean(y_tilde(:, idx(1:MU)), 2);
        z_rec = mean(z_tilde(:, idx(1:MU)), 2);
        
        %% Update sigma
        sigma_old = sigma;
        C = C_handle(N,MU);
        D = D_handle(N,MU);
        s = (1-C)*s + sqrt(MU*C*(2-C))*z_rec;
        
        if CSA_MODE==0
            sigma = sigma * exp((norm(s)^2-E_chi) / (2*D*E_chi)); %HGB, AR
        elseif CSA_MODE == 11 || CSA_MODE == 12
            sigma = sigma * exp(1/D * (norm(s)/E_chi - 1)); % Hansen deradomized
        elseif CSA_MODE == 2
            sigma = sigma * exp(C/D * (norm(s)/E_chi - 1)); % Hansen tutorial
        end     

        %% g>=1: methods write at generation g+1
        if strcmpi(pc_input.pc_method, 'PSA')
            [PSA, p_t, delta_theta, norm2_p_dm, norm2_p_dS] = ...
                pc_psa_simpli(g, y_rec, y_old, sigma, sigma_old, p_t, delta_theta, norm2_p_dm, norm2_p_dS, pc_input, C, D, E_chi, MU); 
        elseif strcmpi(pc_input.pc_method, 'ROSAES')
            [dF, ROSAES] = pc_rosaes(pc_input, g, dF, F_mean_g);
        elseif strcmpi(pc_input.pc_method, 'QGAIN')
            [Q, QGAIN] = pc_qgain(pc_input, g, Q, F_med_g, f_tilde, LAM, N);
        elseif strcmpi(pc_input.pc_method, 'LENGTHZ')
            [Z, LENGTHZ, P_LENGTHZ(g)] = pc_lengthz(pc_input, g, Z, z_rec, MU, N);
        else
            PSA = 0;
            LENGTHZ = 0;
            ROSAES = 0;
            QGAIN = 0;      
        end

        %% g>=L: methods write at generation g+1
        if g>=pc_input.L  % number of fitness values
            debug_g = nan; % nan to deactivate plot
            if strcmpi(pc_input.pc_method, 'PCCMSA')
                [PH, H, PCCMSA] = pc_pccmsa(pc_input, g, PH, H, F_rec_g, debug_g); % WHICH F-values?
            elseif strcmpi(pc_input.pc_method, 'RANK')
                [P_RANK, RANK] = pc_rank(pc_input, g, P_RANK, f_g);
            elseif strcmpi(pc_input.pc_method, 'INFNOISE')
                error('check infnoise');
                %[D2, P_INFNOISE, INFNOISE] = pc_infnoise(pc_input, g, D2, P_INFNOISE, z_tilde(:, idx(1:MU)), MU, N);
            elseif strcmpi(pc_input.pc_method, 'HYBRID')
                error('check hybrid');
            end
        else
            PCCMSA = 0;
            RANK = 0;
            INFNOISE = 0;
            HYBRID = 0;
        end

        %% g>=L: APOP evaluates fitness differences
        if g>=pc_input.L
            if strcmpi(pc_input.pc_method, 'APOP')
                if g>=150
                    dummy = 1;
                end
                [N_up, APOP] = pc_apop(pc_input, g, N_up, F_med_g);
                % [N_up, APOP] = pc_apop_decorr(pc_input, g, N_up, F_med_g);
            end
        else
            APOP = 0;
        end
        
        %% Control
        mu_old = MU;
        P = eval(pc_input.pc_method); % Assign Performance depending on chosen PCS
        if pc_input.pc_on == 1 && wait <= 0
            if strcmpi(pc_input.pid_method, 'pid_v1') %increase/decrease method
                if P==0     % neutral performance
                    mu_new = MU;
                elseif P>0  % positive performance
                    mu_new = MU*pc_input.fac_dec;
                    mu_new = floor(mu_new);
                elseif P<0  % negative performance
                    mu_new = MU*pc_input.fac_inc;
                    mu_new = ceil(mu_new);
                end
            elseif strcmpi(pc_input.pid_method, 'pid_v2')
                mu_new = MU*exp(-pc_input.damp*P);
            elseif strcmpi(pc_input.pid_method, 'pid_psa')
                assert(strcmpi(pc_input.pc_method, 'PSA')) %not defined for other methods
                mu_new = MU*exp(pc_input.beta_psa*(1-(norm2_p_dm(g+1)+norm2_p_dS(g+1))/pc_input.alpha_psa ));
            end
            % keep in bounds and round
            mu_new = min(max(mu_new, pc_input.mu_min), pc_input.mu_max);

            % new mu, adjust params
            if mu_new~=MU 
                wait = pc_input.wait; % wait only if pop. change happened
                if pc_input.sigma_rescale ~= 0
                    sigma = sigma*(mu_new/MU)^(pc_input.sigma_rescale);  
                end
                MU = mu_new;
                LAM = floor(MU/THETA);
            end
        else % no change, decrease wait
            wait = wait-1;
        end

        % if pc_input.sigma_rescale ~= 1 && g==1, warning('sigma-change'); end 
         
        %% SAVE DATA using mu_old; value at g+1 obtained using old MU
        F_rec_g(g+1) = FIT.get_f(y_rec, g, norm(y_rec-FIT.y_hat)); %HEMI %fit(y, sigma_eps_g(g), R_noise); median(fit(y_n_mu, sigeps_g(g)));
        F_med_g(g+1) = median(f_tilde(1, idx(1:mu_old))); %APOP
        F_mean_g(g+1) = mean(f_tilde(1, idx(1:mu_old))); %ROSAES
  
        f_g(g+1) = FIT.get_f(y_rec); % noise-removed fitness of recombinant
        sigma_g(g+1) = sigma;
        r_g(g+1) = norm(y_rec-FIT.y_hat);
        %y_g(:,g+1) = y;
        mu_g(g+1) = MU;
        s_g(:,g+1) = s;
        z_g (:,g+1) = z_rec;

        if verbose == 1
            fprintf('Gen: %i, fit: %d, sigma: %d \n', g+1, f_g(g+1), sigma_g(g+1));
        end
        if sigma_g(g+1) < SIGMA_STOP
            fprintf('\t\t SIGMA_STOP. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g+1, f_g(g+1), r_g(g+1), sigma_g(g+1))
            break
        end
        if r_g(g+1) < R_STOP || f_g(g+1) < F_STOP
            fprintf('\t\t F_R_STOP. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g+1, f_g(g+1), r_g(g+1), sigma_g(g+1))
            break
        end
        if feval_budget<=0
            fprintf('\t\t FEVAL_STOP. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g+1, f_g(g+1), r_g(g+1), sigma_g(g+1))
            break           
        end

    end
    if g+1 == G
        fprintf('\t\t N_G limit reached. Gen: %i, fit: %d, delta_r: %d, sigma %d \n', g+1, f_g(g+1), r_g(g+1), sigma_g(g+1));
    end 
   
    %% Truncate
    g_end = g+1;
    f_g = f_g(1:g_end);
    r_g = r_g(1:g_end);
    sigma_g = sigma_g(1:g_end);
    pc_data.mu_g = mu_g(1:g_end);  

    if strcmpi(pc_input.pc_method, 'PCCMSA')
        pc_data.PH = PH(1:g_end);
        pc_data.F_eval = F_rec_g(1:g_end);
    elseif strcmpi(pc_input.pc_method, 'PSA')
        pc_data.norm2_p_dm = norm2_p_dm(1:g_end);
        pc_data.norm2_p_dS = norm2_p_dS(1:g_end);
    elseif strcmpi(pc_input.pc_method, 'APOP')
        pc_data.N_up = N_up(1:g_end);
    elseif strcmpi(pc_input.pc_method, 'INFNOISE')
        pc_data.D2 = D2(1:g_end);
        pc_data.P_INFNOISE = P_INFNOISE;
    elseif strcmpi(pc_input.pc_method, 'ROSAES')
        pc_data.dF = dF(1:g_end);
    elseif strcmpi(pc_input.pc_method, 'QGAIN')
        pc_data.Q = Q(1:g_end);
    elseif strcmpi(pc_input.pc_method, 'RANK')
        pc_data.P_RANK = P_RANK(1:g_end);
    elseif strcmpi(pc_input.pc_method, 'LENGTHZ')
        pc_data.Z = Z(1:g_end);
        pc_data.P_LENGTHZ = P_LENGTHZ(1:g_end);
    elseif strcmpi(pc_input.pc_method, 'HYBRID')
        pc_data.P_RANK = P_RANK(1:g_end);
        pc_data.P_INFNOISE = P_INFNOISE(1:g_end);
    end

    %% PSA dot-produt of consecutive steps
    % for i=1:g_end-1
    %     sp(i) = p_t(1:N,i)'*z_g(:,i+1);
    % end
    % figure; hold on; 
    %     plot(sp)
    %     yline(mean(sp), 'r--', 'Alpha',1, 'LineWidth',1)
    % dummy=1;

end


    