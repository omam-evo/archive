
clear

%% CSA 1
data = {'const_pcOFF_csa11','const_pcOFF_csa12','const_pcOFF_csa2',...
    'adapt_csa11_v2_APOP','adapt_csa11_v2_PCCMSA','adapt_csa2_v2_PSA'};

%% CSA1 vs CSA2 const. pop.
% data = {'tr50_pcOFF_csa1', 'tr50_pcOFF_csa2'};

%% SWEEP CSA1
% data = dir('tr50_csa1_*');
% data = {'tr50_pcOFF_csa1',data.name};

%% SWEEP CSA2
% data = dir('tr50_csa2_*');
% data = {'tr50_pcOFF_csa2',data.name};

% color = [0,0,0; 0,0,1; 0,1,0; 1,0,0; 0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980];
% style = {'o-','o-','o-','o-','x-.','x-.','x-.'};
color = [1,0,0; 0,0.9,0; 0,0,1; 1,0.5,0; 0 1 1; 1 0 1; 0.8500 0.3250 0.0980];
style = {'o--','*-.','x:','square-','diamond-','^-','x:'};

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;

for i=1:length(data)
    load(fullfile(data{i},'data.mat'), 'pc_summary');

    myplot(f1,f2,f3,f4, pc_summary, data{i}, color(i,:), style{i});
    
    % if i==1
    %     ert = figure;
    % else
    %     figure(ert);
    % end    
    % myplot("$E_r$", pc_summary, data{i}, color(i,:), style{i});

    clear pc_summary
end
    
w = 6;
h = 4;
s = 8;

figure(f1); %PS
    set(gca, 'XScale','log'); 
    xlabel('Variation of $N$'); yticks([0:0.2:1]);
    myfigsize(gcf,w,h,s,s);
    exportgraphics(gcf,'ps_var1.pdf', 'ContentType','vector');
    savefig(gcf,'ps_var1.fig');
figure(f2);  %ERT
    set(gca, 'XScale','log','YScale','log');
    xlabel('Variation of $N$');
    ylim([1e4,1e8]); yticks([10.^(1:1:20)]);
    myfigsize(gcf,w,h,s,s);
    exportgraphics(gcf,'er_var1.pdf', 'ContentType','vector');
    savefig(gcf,'er_var1.fig');
figure(f3);  %PS VAR2
    set(gca, 'XScale','log');
    xlabel('Variation of $N$, $A$'); yticks([0:0.2:1]);
    myfigsize(gcf,w,h,s,s);
    exportgraphics(gcf,'ps_var2.pdf', 'ContentType','vector');
    savefig(gcf,'ps_var2.fig');
figure(f4); %ERT VAR2
    set(gca, 'XScale','log','YScale','log');
    xlabel('Variation of $N$, $A$');
    ylim([1e4,1e8]); yticks([10.^(1:1:20)]);
    myfigsize(gcf,w,h,s,s);
    exportgraphics(gcf,'er_var2.pdf', 'ContentType','vector');
    savefig(gcf,'er_var2.fig');

function [] = myplot(f1,f2,f3,f4, pc_summary, name, color, style)

    id_repeat = 5;
    N_LIST = [10,30,100,300,1000];
    A_LIST = [65,33,12,7,3];
    ALPHA_LIST = [10*pi,7*pi,4*pi,3*pi,2*pi];
    markersize = 6;

    rep_ps = pc_summary.('$P_S$')(id_repeat);
    rep_er = pc_summary.('$E_r$')(id_repeat);

    %% N varied, A=3, ALPHA=2pi
    figure(f1); hold on;
    plot(N_LIST, pc_summary.('$P_S$')(1:5), style, 'MarkerSize', markersize, 'DisplayName', name, 'Color', color);

    figure(f2); hold on;
    plot(N_LIST, pc_summary.('$E_r$')(1:5), style, 'MarkerSize', markersize, 'DisplayName', name, 'Color', color);

    figure(f3); hold on;
    plot(N_LIST, [pc_summary.('$P_S$')(6:9); rep_ps], style,'MarkerSize', markersize, 'DisplayName', name, 'Color', color);

    figure(f4); hold on;
    plot(N_LIST, [pc_summary.('$E_r$')(6:9); rep_er], style,'MarkerSize', markersize, 'DisplayName', name, 'Color', color);


    %% N,A varied, ALPHA=2pi
    % subplot(2,1,2); hold on;
    % set(gca, 'XScale','log');
    % xlabel('Variation of $N$, $A$');
    % 
    % plot(N_LIST, [pc_summary.(field)(6:9); rep], style,'MarkerSize', markersize, 'DisplayName', name, 'Color', color);

    %% N,ALPHA varied, A=3
    % subplot(3,1,3); hold on;
    % set(gca, 'XScale','log');
    % xlabel('Variation of $N$, $\alpha$');
    % 
    % plot(N_LIST, [pc_summary.(field)(10:13); rep], style,'MarkerSize', markersize, 'DisplayName', name, 'Color', color);
end
  
