clear

%% Switch between CSAs via LIST_CSA:
% 'HanV1': (4a),'HanV1_linN': (4b),'HanV2': (4c)
LIST_CSA = {'HanV1'}; %{'HanV1', 'HanV1_linN','HanV2'};  {'SA'}
%% Switch between configurations LIST_CONFIG = [1,2,3,4]
% 1: top left, 2: top right, 3: bottom left, 4: bottom right
LIST_CONFIG = [1]; %[1,2,3,4];
%% Misc:
% > FAC_TRIAL>1 increases number of trials for smoother data
% > Some RES_SIGN_ZERO are given, others are nan; the given values were calculated beforehand with 1-gen experiments
% to improve accuracy; if value is 'nan', sigma^*_0 is calculated from numeric solution of N-dependent progress rate
% figures are saved as '*.fig' and must be merged

MODE = 'CSA';
THETA = 0.5;
FAC_TRIAL = 1;
SAVE_LOAD = 'SAVE';
COLOR = [0,0,1; 1,0,0; 0.4,0.8,0]; 

% G_MAX = @(N) 10000*(N<=100)+2*N*(N>100);
G_MAX = @(N) 100000;
SET_META_EP = 1;

for j = 1:length(LIST_CONFIG)
    CONFIG = LIST_CONFIG(j);
    
    if strcmpi(MODE,'CSA')
        %f1 = figure;
        %f2 = figure;
    end

    if CONFIG == 0 %% N VARIATION
        N_LIST = [100];
        MU_LIST = 10000*ones(1,length(N_LIST));
        TRIAL_LIST = [1]; %1*ones(1,length(N_LIST));
        MODE_MU_N = 1;
        configstr = 'test';

    elseif CONFIG == 1 %% MU var, N small
        MU_LIST = [10,30,100,300,1000,3000,10000];
        N_LIST = 10*ones(1,length(MU_LIST));
        RES_SIGN_OPT = nan*ones(1,length(N_LIST)); %[4.0006    5.4007    8.2508   11.4008   15.6008   18.0009];
        %RES_SIGN_ZERO = [8.9401   16.3441   30.4921   53.0671   97.2401  168.6601 307.3343];
        RES_SIGN_ZERO = [8.87183161936282	16.2921692534914	30.3655289251120	53.1546848900581	96.7832840121050	168.427075505695	306.886686927171];
        TRIAL_LIST = [200,200,200,100,100,100,100]*FAC_TRIAL; %[100,100,100,100,40,40,40];
        MODE_MU_N = 1;
        configstr = 'conf1';

    elseif CONFIG == 2
        MU_LIST = [10,30,100,300,1000,3000,10000];
        N_LIST = 100*ones(1,length(MU_LIST));
        RES_SIGN_OPT = nan*ones(1,length(MU_LIST)); %[6.30055	11.2006	18.55065 27.9007 42.330751 44.25085];
        RES_SIGN_ZERO = nan*ones(1,length(MU_LIST)); %[12.16613	24.92011	47.9651	84.35109	154.70009	268.45009];     
        TRIAL_LIST = [10,10,5,5,3,3,2]*FAC_TRIAL; %[100,100,100,100,40,40,40];
        MODE_MU_N = 1;
        configstr = 'conf2';

    elseif CONFIG == 3 %% N VARIATION
        N_LIST = [10,20,50,100,200,500,1000]; %,500,1000
        MU_LIST = 1000*ones(1,length(N_LIST));
        RES_SIGN_OPT = nan*ones(1,length(N_LIST)); %[15.6008500000000	17.8508500000000	29.2008000000000	34.0008000000000	50.2507500000000	62.5007500000000	88.5007000000000];
        RES_SIGN_ZERO = [96.8241	109.5991	132.5681	nan nan nan nan];
        TRIAL_LIST = [10,10,10,5,2,1,1]*FAC_TRIAL; %[100,100,100,50,25,10,5]; %1*ones(1,length(N_LIST));
        MODE_MU_N = 2;
        configstr = 'conf3';

    elseif CONFIG == 4
        N_LIST = [10,20,50,100,200,500,1000];
        MU_LIST = 2*N_LIST;
        RES_SIGN_OPT = nan*ones(1,length(MU_LIST));
        RES_SIGN_ZERO = [13.122	21.393	41.295 nan nan nan nan];
        TRIAL_LIST = ceil([10,10,10,5,2,1,0.5]*FAC_TRIAL); %[100,100,100,50,25,10,5];
        MODE_MU_N = 2;
        configstr = 'conf4';

    elseif CONFIG == 5
        MU_LIST = [10,30,100,300,1000,3000];
        N_LIST = 1000*ones(1,length(MU_LIST));
        RES_SIGN_OPT = nan*ones(1,length(MU_LIST)); %[7.65055	18.45055	35.20060	55.65065	88.50070	128.25075];
        RES_SIGN_ZERO = nan*ones(1,length(MU_LIST)); %[14.705135	36.9821	79.3761	144.37209	267.86009	466.83009];
        TRIAL_LIST = 10*ones(1,length(N_LIST));
        MODE_MU_N = 1;
        configstr = 'conf5';
    end

    for i = 1:length(LIST_CSA)
        CSA_TYPE = LIST_CSA{i};
        savestr = [CSA_TYPE,'_conf',num2str(CONFIG),'.mat'];

        if strcmpi(SAVE_LOAD,'SAVE')
            [sign_ss_res, ss_path, C_handle, D_handle, r_g, sigma_g] = fct_get_ss(THETA, MODE,CSA_TYPE,N_LIST,MU_LIST,TRIAL_LIST,G_MAX,SET_META_EP);
            save(savestr);
        else
            load(fullfile('data',savestr),'sign_ss_res','C_handle','D_handle');
        end     
        
        SET_B_SIMPLI = 0;
        cv = e_vartheta_a_b(THETA,1,0);
        sph = Sphere;
        sign_pred = nan*MU_LIST;
        NUM_SIGN_OPT = nan*MU_LIST;
        NUM_SIGN_ZERO = nan*MU_LIST;
        gamma = nan*MU_LIST;

         
        for m=1:length(MU_LIST)     
            MU = MU_LIST(m);
            N = N_LIST(m);
            LAM = round(MU/THETA); 
            C_MU_LAM = e_mu_lam_a_b_v2(MU,LAM,1,0);
            if strcmpi(MODE,'CSA')
                if strcmpi(CSA_TYPE,'HanV2')
                    C = C_handle(N,MU);
                    D = 1 + 1/C_handle(N,MU) + 2*max(0, sqrt((MU-1)/(N+1))-1)/C_handle(N,MU);
                else %HanV1
                    C = C_handle(N,MU);
                    D = D_handle(N,MU);
                end

                sign_pred(m) = csa_predict_sign(MU,LAM,N,C,D,SET_B_SIMPLI);
                gamma(m) = csa_get_gamma(THETA,C,D,N);
            end
            [NUM_SIGN_OPT(m),NUM_SIGN_ZERO(m)] = sph.signorm_opt_sph(C_MU_LAM, MU, LAM, N, {'full'});
            
        
        end
        if MODE_MU_N == 1
            x = MU_LIST;
            xstr = '$\mu$';
        else
            x = N_LIST;
            xstr = '$N$';
        end

        for k=1:length(N_LIST)
            if isnan(RES_SIGN_ZERO(k))
                SIGN_ZERO(k) = NUM_SIGN_ZERO(k);
            else
                SIGN_ZERO(k) = RES_SIGN_ZERO(k);
            end
            if isnan(RES_SIGN_OPT(k))
                SIGN_OPT(k) = NUM_SIGN_OPT(k);
            else
                SIGN_OPT(k) = RES_SIGN_OPT(k);
            end
        end
        if length(MU_LIST)>1
            if strcmpi(MODE,'CSA')
                % figure(f1);
                %     hold on;
                %     plot(x, sign_ss_res, '-', 'color', COLOR(i,:), 'DisplayName', CSA_TYPE)
                %     plot(x, sign_pred, '--', 'color', COLOR(i,:), 'DisplayName', CSA_TYPE)
                % 
                %     if i==length(LIST_CSA) % last config of CSA
                %         plot(x, SIGN_ZERO, 'k-.');
                %         xlabel(xstr);
                %         ylabel('$\sigma^*_\mathrm{ss}$')
                %         xlim([x(1),x(end)]);
                %         set(gca,'Xscale','log')
                %         set(gca,'Yscale','log')
                %     end
    
                %% Ratio
                if i==1
                    f2 = figure; hold on;
                end
                    sign_ratio = sign_ss_res./SIGN_ZERO;
                    plot(x, sign_ratio, '-', 'color',COLOR(i,:), 'DisplayName', CSA_TYPE);
                    plot(x, gamma, '--', 'color', COLOR(i,:), 'DisplayName', CSA_TYPE)
                if i==length(LIST_CSA)
                    plot(x, SIGN_OPT./SIGN_ZERO, 'k:'); 
                    xlabel(xstr);
                    xlim([x(1),x(end)]);
                    ylabel('$\gamma$');
                    set(gca,'Xscale','log');
                end
            else % SA (not CSA)
                figure; hold on;
                    sign_ratio = sign_ss_res./SIGN_ZERO;
                    plot(x, sign_ratio, '-', 'color',COLOR(i,:), 'DisplayName', MODE);
                    if i==length(LIST_CSA)
                        plot(x, SIGN_OPT./SIGN_ZERO, 'k:'); 
                        xlabel(xstr);
                        xlim([x(1),x(end)]);
                        ylabel('$\gamma$');
                        set(gca,'Xscale','log');
                    end
            end
                
        end
    end % CSA

    %% Progress rate of each config
    % colors_phi = lines(length(MU_LIST));
    % f3 = figure; hold on; legend;
    % for m=1:length(MU_LIST)
    %     MU = MU_LIST(m);
    %     N = N_LIST(m);
    %     name = ['$\mu$=',num2str(MU), ' $N$=', num2str(N)];
    % 
    %     SIGNZERO = (8*N)^(1/4)*sqrt(cv*MU);
    %     sph = Sphere;
    %     x = linspace(0,SIGNZERO*1.1,101);
    %     phi = sph.phi(x, cv, MU, LAM, N, {'full'});
    % 
    %     legend('AutoUpdate','on')
    %     plot(x, phi, '-', 'color', colors_phi(m,:), 'DisplayName',name);
    %     legend('AutoUpdate','off')
    %     plot(x, cv*sqrt(2*N)-x.^2/(2*MU), '--', 'color', colors_phi(m,:));
    %     xlabel('$\sigma^*$');
    %     ylabel('$\varphi^*$');
    %     %xline(sign_ss_res(m), '-')
    %     %xline(sign_pred(m), '--')
    % end

    if strcmpi(MODE,'CSA')
        fs = 9;
        %myfigsize(f1,8,5,fs,fs);

        figure(f2);
        ylim([0.8,1]);
        
        myfigsize(f2,8,5,fs,fs);
        
        % myfigsize(f3,8,5,fs,fs);
        % set(gca,'Xscale','log')
        % set(gca,'Yscale','log')
        
        %savefig(f1,['data/sign_',configstr, '.fig']); 
        %exportgraphics(f1,['data/sign_',configstr, '.pdf'],'ContentType','vector');
        savefig(f2,['gamma_', configstr, '.fig']); 
        %exportgraphics(f2,['gamma_', configstr, '.pdf'],'ContentType','vector');
        % savefig(f3,['data/phi_', configstr, '.fig']); 
        % figure(f3); legend('off');
        % exportgraphics(f3,['data/phi_', configstr, '.pdf'],'ContentType','vector');
    end

end %config


