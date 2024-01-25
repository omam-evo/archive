function [] = get_scaling_lines(N_LIST, base)
    hold on;
    facs = [25,50,75,100,150,300];
    x = min(N_LIST):0.01:max(N_LIST);
    for p=1:length(facs)
        fac = facs(p);
        if p==1
            legend('autoupdate','on');
        else
            legend('autoupdate','off');
        end
        plot(x, fac*x, 'c:', 'DisplayName', '$N$');
    end
end