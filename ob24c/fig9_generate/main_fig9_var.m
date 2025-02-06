clear
addpath('pc_input');

%% To generate Tables II, III, and IV:
% > execute script 

%% Describe variation; constant params
prefix = 'ras';
LIST_METHOD = {'APOP'}; %{'APOP', 'PSA', 'PCCMSA'};
LIST_CSA = [11]; %[0,11,12,2];  %[0: Arnold, 11: HanSqrtN, 12: HanLinN, 2: HanNew]
LIST_CONFIG = [2];

%% OFF
% pc_input.mu_max = 10;
% pc_input.mu_min = 10;
% pc_input.pc_on = 0; prefix = [prefix, '_pcOFF'];
%% ON
pc_input.mu_max = 1024;
pc_input.mu_min = 4;
pc_input.pc_on = 1;

%% Additional Params
pc_input.TRIALS = 50;
pc_input.SEED = 1; warning('seed');

pc_input.theta = 0.5;
pc_input.VAR_TYPE = 1;
MU = pc_input.mu_min;
LAM = pc_input.mu_min/pc_input.theta;


%% VARIATION 
% [fconf, init_stop, readme] = get_sphere_var(MU, LAM);
% [fconf, init_stop, readme] = get_random_var();
% [fconf, init_stop, readme] = get_sphere_hgbTest(100);

% [fconf, init_stop, readme] = get_rastrigin_varN(MU, LAM);
% [fconf, init_stop, readme] = get_rastrigin_varA(MU, LAM, 0);
% [fconf, init_stop, readme] = get_rastrigin_varAlpha(MU, LAM, 0);

% ----- Sweep all Rastrigin
[fconf, init_stop, readme] = get_rastrigin_varN(MU, LAM);
[fconf(6:9), init_stop(6:9), readme(6:9)] = get_rastrigin_varA(MU, LAM, 5); % remove 5th config
% [fconf(10:13), init_stop(10:13), readme(10:13)] = get_rastrigin_varAlpha(MU, LAM, 5); % remove 5th config

% SLOW ADAPTATION USING CSA12
% ----- Sweep Sphere and Random
% [fconf, init_stop, readme] = get_sphere_var_CSA12(MU, LAM);
% [fconf(4:6), init_stop(4:6), readme(4:6)] = get_random_var_CSA12();

%% MAIN
tic
cc = 1;
for s=1:length(LIST_CSA)
    fprintf('----- CSA=%i -----\n', LIST_CSA(s));
    pc_input.CSA_MODE = LIST_CSA(s); prefix_s = [prefix, '_csa',num2str(pc_input.CSA_MODE)];

    for c=1:length(LIST_CONFIG)

        config = LIST_CONFIG(c);
        fprintf('\t ----- config=%i -----\n', config);

        switch config
            case 1
                pc_get_config = @(pc_input, N) input_v1(pc_input, N); prefix_c = [prefix_s, '_v1'];
            case 2
                pc_get_config = @(pc_input, N) input_v2(pc_input, N); prefix_c = [prefix_s, '_v2'];
            case 3
                pc_get_config = @(pc_input, N) input_v3(pc_input, N); prefix_c = [prefix_s, '_v3'];
            case 4
                pc_get_config = @(pc_input, N) input_v4(pc_input, N); prefix_c = [prefix_s, '_v4'];
        end

        for m = 1:length(LIST_METHOD)

            pc_input.pc_method = LIST_METHOD{m}; descr = [prefix_c, '_', pc_input.pc_method];  %'APOP', 'PSA', 'PCCMSA'
            
            [data_input{cc}, data_summary{cc}] = main_pc(pc_input, fconf, init_stop, readme, pc_get_config, descr);

            cc = cc+1;
        end % method
    end % input config
end %CSA variant
toc
        

