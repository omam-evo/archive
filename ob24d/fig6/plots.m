addpath('../code')

ERT_log = log10(ERT);
[maxval, max_id] = max(ERT_log, [], "all", "linear");
[minval, min_id] = min(ERT_log, [], "all", "linear");

[~,t_max,v_max,m_max] = ind2sub(size(ERT), max_id);
[~,t_min,v_min,m_min] = ind2sub(size(ERT), min_id);

global SAVE w h f plot_minval plot_maxval ert_levels ert_colors caxis_ticks lam_ticks prob_ticks ps_levels ps_nvals ps_colors

%% GENERATE AND SAVE
SAVE = 1;
ids_theta = [3,5] %[3,5];
ids_tau = [1,4] %[1,4];
ids_lam = [5];

%% ERT SETTINGS N50
assert(N==50);
w = 6;
h = 5;
f = 8;
loc = 'eastoutside';
plot_minval = 5.25;
plot_maxval = 8;
ert_levels = [plot_minval:0.25:plot_maxval];
caxis_ticks = ert_levels;
ert_colors = length(ert_levels)-1;
lam_ticks = [200,500,1000,1500,2000];

%% ERT SETTINGS N200
% assert(N==200);
% w = 6;
% h = 5;
% f = 8;
% loc = 'eastoutside';
% plot_minval = 4.75;
% plot_maxval = 7;
% ert_levels = [plot_minval:0.25:plot_maxval];
% caxis_ticks = [5:0.5:7];
% ert_colors = length(ert_levels)-1;
% lam_ticks = [200,400,600,800,1000];

%% PS SETTINGS N50
prob_ticks = 0:0.2:1;
% ps_levels = [0.0:0.1:0.9,0.95];
ps_levels = [0:0.1:0.9,1.05];
ps_nvals = length(ps_levels)-1;
ps_colors = hot(ps_nvals);
ps_colors(2:7,:) = ps_colors(1:6,:);
ps_colors(1,:) = [0,0,0];

%% LAM VS. TAU
for i=1:length(ids_theta)

    id_theta = ids_theta(i); 
    savestr = ['vt',strrep(num2str(THETA_LIST(id_theta)),'.','p')];
    titlestr = ['$\vartheta$=',num2str(THETA_LIST(id_theta))];

    [x,y] = meshgrid(LAM_LIST, TAU_LIST);
    z = 0*x;
    p = 0*x;
    for t=1:N_tau
        for v=id_theta %1:N_theta
            for m=1:N_lam
                z(t,m) = ERT_log(1,t,v,m);
                p(t,m) = P_S(1,t,v,m);
            end
        end
    end

    ert_plot(x,y,z, '$\lambda$', '$\tau$', lam_ticks, TAU_LIST, titlestr, savestr);
    ps_plot(x,y,p, '$\lambda$', '$\tau$', lam_ticks, TAU_LIST, titlestr, savestr);

    % figure; hold on;    
    %     xlabel('$\lambda$'); ylabel('$\tau$')
    %     contourf(x,y,z,'LevelList',ert_levels);
    %     xticks(lam_ticks);
    %     yticks(TAU_LIST);
    %     cc = parula(ert_colors);
    %     cc = flipud(cc);
    %     colormap(cc);
    %     clim([plot_minval, plot_maxval]);
    %     if SAVE==0
    %         title(['$\vartheta$=',num2str(THETA_LIST(id_theta))]);
    %         cb = colorbar('XTick', caxis_ticks,'location',loc);          
    %     end
    % myfigstyle(gcf,w,h,f,f);
    % if SAVE==1
    %     savefig([savestr,'_ert.fig']);
    %     exportgraphics(gca,[savestr,'_ert.pdf']);
    % end
    % 
    % figure; hold on;
    % 
    %     xlabel('$\lambda$'); ylabel('$\tau$')
    %     contourf(x,y,p,'LevelList',ps_levels);
    %     xticks(lam_ticks);
    %     yticks(TAU_LIST);
    %     colormap(ps_colors);
    %     clim([0, 1]);
    %     if SAVE==0
    %         title(['$\vartheta$=',num2str(THETA_LIST(id_theta))]);   
    %         colorbar('XTick', prob_ticks,'location',loc);
    %     end   
    % myfigstyle(gcf,w,h,f,f);
    % if SAVE==1
    %     savefig([savestr,'_ps.fig']);
    %     exportgraphics(gcf,[savestr,'_ps.pdf']);
    % end
end

%% MU VS. THETA
for i=1:length(ids_tau)

    id_tau = ids_tau(i); 
    savestr = ['tau',strrep(num2str(TAU_LIST(id_tau)),'.','p')];
    titlestr = ['$\tau$=',num2str(TAU_LIST(id_tau))];
    
    [x,y] = meshgrid(LAM_LIST, THETA_LIST);
    z = 0*x;
    p = 0*x;
    for t=id_tau %1:N_tau
        for v=1:N_theta
            for m=1:N_lam
                z(v,m) = ERT_log(1,t,v,m);
                p(v,m) = P_S(1,t,v,m);
            end
        end
    end
    ert_plot(x,y,z, '$\lambda$', '$\vartheta$', lam_ticks, THETA_LIST, titlestr, savestr);
    ps_plot(x,y,p, '$\lambda$', '$\vartheta$', lam_ticks, THETA_LIST, titlestr, savestr);

    % figure; hold on;       
    %     xlabel('$\lambda$'); ylabel('$\vartheta$');
    %     contourf(x,y,z,'LevelList',ert_levels);
    %     xticks(lam_ticks);
    %     yticks(THETA_LIST);
    %     cc = parula(ert_colors);
    %     cc = flipud(cc);
    %     colormap(cc);
    %     clim([plot_minval, plot_maxval]);
    %     if SAVE==0
    %         title(['$\tau$=',num2str(TAU_LIST(id_tau))]);
    %         cb = colorbar('XTick', caxis_ticks,'location',loc); 
    %     end
    % myfigstyle(gcf,w,h,f,f);
    % if SAVE==1
    %     savefig([savestr,'_ert.fig']);
    %     exportgraphics(gca,[savestr,'_ert.pdf']);
    % end
    % 
    % figure; hold on;       
    %     xlabel('$\lambda$'); ylabel('$\vartheta$');
    %     contourf(x,y,p,'LevelList',ps_levels);
    %     xticks(lam_ticks);
    %     yticks(THETA_LIST);
    %     colormap(ps_colors);
    %     clim([0, 1]);
    %     if SAVE==0
    %         title(['$\tau$=',num2str(TAU_LIST(id_tau))]);
    %         colorbar('XTick', prob_ticks,'location',loc);     
    %     end
    % myfigstyle(gcf,w,h,f,f);
    % if SAVE==1
    %     savefig([savestr,'_ps.fig']);
    %     exportgraphics(gcf,[savestr,'_ps.pdf']);
    % end
end


%% THETA VS. TAU
for i=1:length(ids_lam)

    id_lam = ids_lam(i); 
    savestr = ['lam',strrep(num2str(LAM_LIST(id_lam)),'.','p')];
    titlestr = ['$\lambda$=',num2str(LAM_LIST(id_lam))];

    [x,y] = meshgrid(TAU_LIST, THETA_LIST);
    z = 0*x;
    p = 0*x;
    for t=1:N_tau %1:N_tau
        for v=1:N_theta
            for m=id_lam
                z(v,t) = ERT_log(1,t,v,m);
                p(v,t) = P_S(1,t,v,m);
            end
        end
    end
    ert_plot(x,y,z, '$\tau$', '$\vartheta$', TAU_LIST, THETA_LIST, titlestr, savestr);
    ps_plot(x,y,p, '$\tau$', '$\vartheta$', TAU_LIST, THETA_LIST, titlestr, savestr);

    % figure; hold on;        
    %     xlabel('$\tau$'); ylabel('$\vartheta$');
    %     contourf(x,y,z,'LevelList',ert_levels);
    %     xticks(TAU_LIST);
    %     yticks(THETA_LIST);
    %     cc = parula(ert_colors);
    %     cc = flipud(cc);
    %     colormap(cc);   
    %     clim([plot_minval, plot_maxval]);
    %     if SAVE==0
    %         title(['$\lambda$=',num2str(LAM_LIST(id_lam))]);
    %         cb = colorbar('XTick', caxis_ticks,'location',loc); 
    %     end
    % myfigstyle(gcf,w,h,f,f);
    % if SAVE==1
    %     savefig([savestr,'_ert.fig']);
    %     exportgraphics(gca,[savestr,'_ert.pdf']);
    % end
    
    % figure; hold on;     
    %     xlabel('$\tau$'); ylabel('$\vartheta$');
    %     contourf(x,y,p,'LevelList',ps_levels);
    %     xticks(TAU_LIST);
    %     yticks(THETA_LIST);
    %     clim([0, 1]);
    %     colormap(ps_colors);
    %     if SAVE==0
    %         title(['$\lambda$=',num2str(LAM_LIST(id_lam))]);
    %         colorbar('XTick', prob_ticks,'location',loc);
    %     end
    % 
    % myfigstyle(gcf,w,h,f,f);
    % if SAVE==1
    %     savefig([savestr,'_ps.fig']);
    %     exportgraphics(gcf,[savestr,'_ps.pdf']);
    % end
end


%% COLOR BAR
width = 6; %dummy width
figure
    c = colorbar('XTick', ert_levels,'location','SouthOutside', 'AxisLocation','in');
    cc = parula(ert_colors);
    cc = flipud(cc);
    colormap(cc);
    clim([plot_minval, plot_maxval]);
    axis off
    myfigstyle(gcf,width,1.3,8,8);
    pos = c.Position;
    % set(c,'Position',[pos(1) pos(2) pos(3) pos(4)])
    % set(c,'Position',[pos(1) 0.08 pos(3) 0.45])
    set(c,'Position',[0.05 0.15 0.90 0.30])
    if SAVE==1
        savefig(['ert.fig']);
        exportgraphics(gcf,['ert.pdf']);
    end

figure
    ps_levels_show1 = ps_levels;
    ps_levels_show1(end) = 1;
    c = colorbar('XTick', ps_levels_show1,'location','SouthOutside', 'AxisLocation','in');
    colormap(ps_colors);
    clim([0,1]);
    axis off
    myfigstyle(gcf,width,1.3,8,8);
    set(c,'Position',[0.05 0.15 0.90 0.30])
    if SAVE==1
        savefig(['ps.fig']);
        exportgraphics(gcf,['ps.pdf']);
    end

%% ert_plot
function [] = ert_plot(x,y,z, x_lab, y_lab, x_tick, y_tick, titlestr, savestr)
    global SAVE w h f plot_minval plot_maxval ert_levels ert_colors caxis_ticks
    figure; hold on;        
        xlabel(x_lab); ylabel(y_lab);
        contourf(x,y,z,'LevelList',ert_levels);
        xticks(x_tick); yticks(y_tick);
        cc = parula(ert_colors);
        cc = flipud(cc);
        colormap(cc);   
        clim([plot_minval, plot_maxval]);
        if SAVE==0
            title(titlestr);
            cb = colorbar('XTick', caxis_ticks); 
        end
    myfigstyle(gcf,w,h,f,f);
    if SAVE==1
        savefig([savestr,'_ert.fig']);
        exportgraphics(gca,[savestr,'_ert.pdf']);
    end
end

%% ps_plot
function [] = ps_plot(x,y,p, x_lab, y_lab, x_tick, y_tick, titlestr, savestr)
global SAVE w h f prob_ticks ps_levels ps_nvals ps_colors
    figure; hold on;     
        xlabel(x_lab); ylabel(y_lab);
        contourf(x,y,p,'LevelList',ps_levels);
        xticks(x_tick); yticks(y_tick);
        clim([0, 1]);
        colormap(ps_colors);
        if SAVE==0
            title(titlestr);
            colorbar('XTick', prob_ticks);
        end
        
    myfigstyle(gcf,w,h,f,f);
    if SAVE==1
        savefig([savestr,'_ps.fig']);
        exportgraphics(gcf,[savestr,'_ps.pdf']);
    end
end

