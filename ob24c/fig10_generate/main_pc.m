%% CONFIG: no var
% RUN_VAR = 0;

%% CONFIG: var
function [pc_input, pc_summary] = main_pc(pc_input, fconf, init_stop, readme, pc_get_config, descr)
RUN_VAR = 1;

if RUN_VAR == 0
    clear;
    addpath('../n0_generic');
    addpath('pc_input');
    addpath('C:/Users/omam/optimize/mat_v2/fit');

    %% OFF theta=1/2
    % pc_input.mu_max = 1024;
    % pc_input.mu_min = 1024;
    % pc_input.pc_on = 0;
    %% ON theta=1/2
    pc_input.mu_max = 1024;
    pc_input.mu_min = 4;
    pc_input.pc_on = 1;

    pc_input.TRIALS = 20;
    pc_input.CSA_MODE = 1;
    pc_input.pc_method = 'APOP';
    pc_input.theta = 0.5;
    pc_input.VAR_TYPE = 1;

    %% PC Config
    pc_get_config = @(pc_input, N) input_slow(pc_input, N);
    % pc_get_config = @(pc_input, N) input_slow_rescale(pc_input, N);
    % pc_get_config = @(pc_input, N) input_fast_wait(pc_input, N);
    
    %% ----- VARIATION -----
    MU = pc_input.mu_min;
    LAM = pc_input.mu_min/pc_input.theta;

    % [fconf, init_stop, readme] = get_sphere_var(MU, LAM);
    % [fconf, init_stop, readme] = get_random_var();
    [fconf, init_stop, readme] = get_rastrigin_varA(MU, LAM);
    % [fconf, init_stop, readme] = get_rastrigin_varN(MU, LAM);
    % [fconf, init_stop, readme] = get_rastrigin_varAlpha(MU, LAM);
    % --- choose configs ---
    id = 1; fconf = fconf(id); init_stop = init_stop(id); readme = readme(id);
    
    % ----- SINGLE EVALUATION -----
    % N = 100;input
    % A = 12; ALPHA = 2*pi;
    % [fconf, init_stop, readme] = get_rastrigin(N, A, ALPHA, MU, LAM);
    % % [fconf, init_stop, readme] = get_random(N);

end

%% SET FROM INPUT
MU = pc_input.mu_min;
LAM = pc_input.mu_min/pc_input.theta;
M = pc_input.TRIALS;

%% GENERAL
PLOT = 1;

%% METHOD VARIATION => NOT YET IMLPEMENTED
% pc_input.VAR_TYPE = 2; 
% VAR_METHOD = 'PC';
% VAR_PC = {'PCCMSA','APOP'}; %'APOP', 'PSA', 'PCCMSA', 'LENGTHZ' 'RANK', 'HYBRID', 'INFNOISE', 'QGAIN'
% NUM_VAR = length(VAR_PC);

%% VARIATION INIT
NUM_VAR = length(fconf);
FEVAL_v_t = nan*zeros(NUM_VAR,M); 
MUMED_v_t = nan*zeros(NUM_VAR,M); 
SUCC_v_t = zeros(NUM_VAR,M);
MU_v = zeros(NUM_VAR,1); 
MU_25 = zeros(NUM_VAR,1); 
MU_75 = zeros(NUM_VAR,1); 
P_S_v = nan*zeros(NUM_VAR,1);
ERT_v = nan*zeros(NUM_VAR,1); 
variation_options = [1,NUM_VAR,pc_input.VAR_TYPE]; %1=j in var-loop

tic
 for j=1:NUM_VAR
    %% VARIATION SETTINGS
    if pc_input.VAR_TYPE == 1
        fit = fconf(j); 
    elseif pc_input.VAR_TYPE == 2
        %pc_input.pc_method = VAR_PC{j}; 
        error('Method variation not implemented.');
    end
    fprintf('pc=%s, fit=%s, N=%i \n', pc_input.pc_method, fit.name, fit.N);

    %% N-dependent quantities
    N = fit.N;
    if pc_input.CSA_MODE == 11 %% CSA Hansen sqrt(N)
        C = @(N,MU) 1/sqrt(N);
        D = @(N,MU) 1/C(N,MU); 
    elseif pc_input.CSA_MODE == 12 %CSA Hansen lin N
        C = @(N,MU) 1/N;
        D = @(N,MU) 1/C(N,MU);        
    elseif pc_input.CSA_MODE == 2 %CSA Hansen new
        C = @(N,MU) (MU+2)/(MU+N+5);
        D = @(N,MU) 1 + C(N,MU) + 2 * max(0, sqrt((MU-1)/(N+1))-1);
    elseif pc_input.CSA_MODE == 0 %HGB, AR
        C = @(N,MU) 1/sqrt(N);
        D = @(N,MU) 1/C(N,MU);        
    end

    if pc_input.CSA_MODE == 0 % HGB, AR
        E_chi = N;
    else % Hansen
        E_chi = sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2)); %sqrt(2)*gamma((N+1)/2)/gamma(N/2)
    end 
    S_0 = 1;

    %% Pop. control
    pc_input = pc_get_config(pc_input, N);

    %% MAIN LOOP
    rng(pc_input.SEED); % evaluation of single config consistent with first eval. over many configs
    for i=1:M
        
        % If needed, implement random y_0 and sigma_0 for each trial

        [y, f_g, r_g, sigma_g, y_g, pc_data] = ...
            csa_es_meas(pc_input, fit, pc_input.CSA_MODE, N, MU, LAM, init_stop(j).Y_0, init_stop(j).SIGMA_0, S_0, C, D, E_chi, init_stop(j).SIGMA_STOP, init_stop(j).R_STOP, init_stop(j).F_STOP, init_stop(j).F_EVAL, init_stop(j).G, 0);

        %% Save data (only last trial preserved)
        dyn_res.r_g = r_g;
        dyn_res.f_g = f_g;
        dyn_res.sigma_g = sigma_g;
        dyn_res.mu_g = pc_data.mu_g;

        if f_g(end) < init_stop(j).F_STOP
            SUCC_v_t(j, i) = 1;
        end
        FEVAL_v_t(j, i) = sum(pc_data.mu_g(1:end-1)*(LAM/MU)); %last feval is not taking place
        MUMED_v_t(j, i) = median(pc_data.mu_g);

        if PLOT == 1
            if i==1 % plot the single run for each var.
                fct_plot(pc_input, dyn_res, pc_data, init_stop, [j,NUM_VAR,pc_input.VAR_TYPE], fit, pc_input.theta);
                if M>1
                    warning('First plot over many trials is shown.')
                end
            end
        end
    end
    N_SUCC = sum(SUCC_v_t(j,:));
    P_S_v(j) = N_SUCC/M;
    ERT_v(j) = sum(FEVAL_v_t(j,:))/N_SUCC;
    MU_v(j) = median(MUMED_v_t(j,:));
    MU_25(j) = prctile(pc_data.mu_g,25); %only last is saved
    MU_75(j) = prctile(pc_data.mu_g,75); %only last is saved

end
toc
% pc_summary = table(num2str(MU_v), num2str(TRIALS*ones(NUM_VAR,1)), num2str(P_S_v,'%.2f'), num2str(sum(FEVAL_v_t,2),'%.2e'), num2str(ERT_v,'%.2e'),...
%     'VariableNames',["MU_MED", "TRIALS","P_S", "FEVAL","ERT"], ...
%     'RowNames', [readme{:}]');


if any(strcmpi({fconf(1,:).name},'Rastrigin'))
    vars = ["$M$","$P_S$", "$\mu_\text{med}$","$F_t$","$E_r$"];
    form = {'%i',1,'%.2f',1,'%i',1, '%.1e',2};
else
    vars = ["$\mu_\text{25}$", "$\mu_\text{med}$","$\mu_\text{75}$","$F_t$"];
    form = {'%i',1,'%i',1,'%i',1,'%.1e',1};
end

if RUN_VAR == 1
    mkdir(descr);

    %% ACTIVATE NUMBER FOR UNIQUE ROW IDENTIFIER (NAMES MUST NOT BE EQUAL)
    for i=1:length(readme)
        readme{i} = ['' ,readme{i}{:}]; %(',num2str(i), ') 
    end

    if any(strcmpi({fconf(1,:).name},'Rastrigin'))
        pc_summary = table(M*ones(NUM_VAR,1), P_S_v, MU_v, sum(FEVAL_v_t,2), ERT_v,...
            'VariableNames',vars, ...
            'RowNames', readme); %[readme{:}]'
    else
        pc_summary = table(MU_25, MU_v, MU_75, sum(FEVAL_v_t,2), 'VariableNames',vars, ...
            'RowNames', readme); %[readme{:}]'
    end
     
    input.data = pc_summary;
    input.dataFormatMode = 'column';
    input.dataFormat = form;
    input.makeCompleteLatexDocument = 0;
    latex = latexTable(input);
    
    fid=fopen(fullfile(descr,'table.asv'),'w');
    [nrows,ncols] = size(latex);
    for row = 1:nrows
        fprintf(fid,'%s\n',latex{row,:});
    end
    fclose(fid);
    fprintf('... table.asv saved in working directory\n');
    
    savefig(gcf, fullfile(descr, 'gcf.fig'));
    save(fullfile(descr, 'data.mat'), 'pc_input', 'pc_summary', 'fconf', 'init_stop', 'readme');
end