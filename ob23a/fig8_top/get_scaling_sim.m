function [] = get_scaling_sim(VAR, MU_givenP, P_target, call, str_var, plot_lin_log)
    sim_line = 'ko-';
    sim_label = 'Sim.';

    if strcmpi(call,'firstCall')
        figure; hold on; 
        xlabel(str_var); 
        ylabel(['$\mu','(P_S\geq',num2str(P_target),')$']);

        legend('autoupdate','on');
        plot(VAR, MU_givenP, sim_line, 'DisplayName', sim_label);
        if strcmpi(plot_lin_log, 'log')
            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'log');
        end
        get_ylim = ylim;
        xlim([VAR(1),VAR(end)]);
        ylim(get_ylim);
    elseif strcmpi(call,'secondCall')
        hold on;
        plot(VAR, MU_givenP, sim_line, 'DisplayName', sim_label);
    end

end