function [] = get_scaling_sim(ALPHA_LIST, MU_givenP, P_target, call)
    sim_line = 'ko-';
    sim_label = 'Sim.';

    if strcmpi(call,'firstCall')
        figure; hold on; 
        xlabel('$\alpha$'); 
        ylabel(['$\mu','(P_S\geq',num2str(P_target),')$']);

        legend('autoupdate','on');
        plot(ALPHA_LIST, MU_givenP, sim_line, 'DisplayName', sim_label);
        set(gca, 'YScale', 'log');
        get_ylim = ylim;
        xlim([ALPHA_LIST(1),ALPHA_LIST(end)]);
        ylim(get_ylim);
        %% Ticks
        xticks([pi/2, pi, 2*pi, 4*pi, 8*pi]);
        xticklabels({'$\pi/2$', '$\pi$', '$2\pi$', '$4\pi$', '$8\pi$'}');
        set(gca,'XMinorTick','off','YMinorTick','on');
    elseif strcmpi(call,'secondCall')
        hold on;
        plot(ALPHA_LIST, MU_givenP, sim_line, 'DisplayName', sim_label);
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
        yticks([1,1e1,1e2,1e3])
    end

end