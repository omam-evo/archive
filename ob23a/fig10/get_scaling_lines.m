function [] = get_scaling_lines(ALPHA_LIST, base)
    hold on;
    facs = [0.001,0.01,0.1,1,10,100];
    x = min(ALPHA_LIST):0.01:max(ALPHA_LIST);
    for p=1:length(facs)
        fac = facs(p);
        if p==1
            legend('autoupdate','on');
        else
            legend('autoupdate','off');
        end
        plot(x, fac*x.^2, 'c:', 'DisplayName', '$N^2$');
    end
end