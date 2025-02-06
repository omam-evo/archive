addpath('../code')

ERT_log = log10(ERT);
[maxval, max_id] = max(ERT_log, [], "all", "linear");
[minval, min_id] = min(ERT_log, [], "all", "linear");

[~,t_max,v_max,m_max] = ind2sub(size(ERT), max_id);
[~,t_min,v_min,m_min] = ind2sub(size(ERT), min_id);

global SAVE w h f plot_minval plot_maxval ert_levels ert_colors caxis_ticks lam_ticks prob_ticks ps_levels ps_nvals ps_colors

%% GENERATE AND SAVE
SAVE = 0;
ids_v1 = [1]; 

%% ERT SETTINGS N50
% assert(N==50);
w = 6;
h = 5;
f = 8;
loc = 'eastoutside';
plot_minval = floor(minval);
plot_maxval = ceil(maxval);
ert_levels = [plot_minval:0.5:plot_maxval];
caxis_ticks = ert_levels;
ert_colors = length(ert_levels)-1;


%% PS SETTINGS N50
prob_ticks = 0:0.2:1;
% ps_levels = [0.0:0.1:0.9,0.95];
ps_levels = [0:0.1:0.9,1.05];
ps_nvals = length(ps_levels)-1;
ps_colors = hot(ps_nvals);
ps_colors(2:7,:) = ps_colors(1:6,:);
ps_colors(1,:) = [0,0,0];

%% 
V1_LIST = N_LIST;
N_v1 = length(V1_LIST);
V2_LIST = A_LIST;
N_v2 = length(V2_LIST);
V3_LIST = ALPHA_LIST;
N_v3 = length(V3_LIST);

for i=1:length(ids_v1)

    v1 = ids_v1(i); 
    savestr = ['v1_',strrep(num2str(V1_LIST(v1)),'.','p')];
    titlestr = ['v1_',num2str(V1_LIST(v1))];

    [x,y] = meshgrid(V2_LIST, V3_LIST);
    z = 0*x;
    p = 0*x;
    succ = 0*x;
    for v2=1:N_v2
        for v3=1:N_v3

                fprintf('A=%i, alpha=%d \n', A_LIST(v2), ALPHA_LIST(v3))
                z(v3,v2) = ERT_log(v1,v2,v3);
                p(v3,v2) = P_S(v1,v2,v3);

                if all([dat_suc(v1,v2,v3,:)]==1)
                    flag = 1;
                elseif any([dat_suc(v1,v2,v3,:)]==-1)
                    flag = -1;
                else
                    flag = 0;
                end

                succ(v3,v2) = flag;
        end
    end

    %ert_plot(x,y,z, '$\alpha$', '$A$', V2_LIST, V3_LIST, titlestr, savestr);
    ps_plot(x,y,p, '$A$', '$\alpha$', V2_LIST, V3_LIST, titlestr, savestr);
    succ_plot(x,y,succ, '$A$', '$\alpha$', V2_LIST, V3_LIST, titlestr, savestr);
end

%% ert_plot
function [] = ert_plot(x,y,z, x_lab, y_lab, x_tick, y_tick, titlestr, savestr)
    global SAVE w h f plot_minval plot_maxval ert_levels ert_colors caxis_ticks
    figure; hold on;        
        xlabel(x_lab); ylabel(y_lab);
        contourf(x,y,z,'LevelList',ert_levels);
        %xticks(x_tick); yticks(y_tick);
        cc = parula(ert_colors);
        cc = flipud(cc);
        colormap(cc);   
        clim([plot_minval, plot_maxval]);
        if SAVE==0
            title(titlestr);
            cb = colorbar('XTick', caxis_ticks); 
        end
        set(gca,"XScale", "Log");
        set(gca,"YScale", "Log");
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
        %xticks(x_tick); yticks(y_tick);
        clim([0, 1]);
        colormap(ps_colors);
        if SAVE==0
            title(titlestr);
            colorbar('XTick', prob_ticks);
        end
        set(gca,"XScale", "Log");
        set(gca,"YScale", "Log");
    myfigstyle(gcf,w,h,f,f);
    if SAVE==1
        savefig([savestr,'_ps.fig']);
        exportgraphics(gcf,[savestr,'_ps.pdf']);
    end
end

%% succ_plot
function [] = succ_plot(x,y,succ, x_lab, y_lab, x_tick, y_tick, titlestr, savestr)
global SAVE w h f
    x = reshape(x,numel(x),1);
    y = reshape(y,numel(y),1);
    succ = reshape(succ,numel(succ),1);
    c = nan*zeros(length(succ),3);
    for i=1:length(succ)
        if succ(i)==1
            c(i,:) = [0,1,0];
        elseif succ(i) == 0
            c(i,:) = [0,0,1];
        elseif succ(i) == -1
            c(i,:) = [1,0,0];
        end
    end

    figure; hold on;     
        xlabel(x_lab); ylabel(y_lab);
        scatter(x,y,[], c);
        if SAVE==0
            title(titlestr);
        end
        set(gca,"XScale", "Log");
        set(gca,"YScale", "Log");
    myfigstyle(gcf,w,h,f,f);
    if SAVE==1
        savefig([savestr,'_succ.fig']);
        exportgraphics(gcf,[savestr,'_succ.pdf']);
    end
end
