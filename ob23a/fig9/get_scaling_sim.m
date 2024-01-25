function [] = get_scaling_sim(A_LIST, MU_givenP, P_target, call)
    sim_line = 'ko-';
    sim_label = 'Sim.';

    if strcmpi(call,'firstCall')
        figure; hold on; 
        xlabel('$A$'); 
        ylabel(['$\mu','(P_S\geq',num2str(P_target),')$']);

        legend('autoupdate','on');
        plot(A_LIST, MU_givenP, sim_line, 'DisplayName', sim_label);
        %set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        get_ylim = ylim;
        xlim([A_LIST(1),A_LIST(end)]);
        ylim(get_ylim);
    elseif strcmpi(call,'secondCall')
        hold on;
        plot(A_LIST, MU_givenP, sim_line, 'DisplayName', sim_label);
    end

end