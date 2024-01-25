function [] = get_scaling_lines(N_LIST, base)
    hold on;
    optical_correction = 1;
    facs = optical_correction*base.^(-2:1:5);
    x = min(N_LIST):0.01:max(N_LIST);
    for p=1:length(facs)
        fac = facs(p);
        if p==1
            legend('autoupdate','on');
        else
            legend('autoupdate','off');
        end
        plot(x, fac*x.^(1/2), 'c:', 'DisplayName', '$\sqrt{N}$');
        %plot(x, fac*sqrt(x).*log(x), 'g-.', 'DisplayName', '$\sqrt{N}\log{N}$');
        plot(x, fac*x.^(5/8), 'g-.', 'DisplayName', '$N$');
        plot(x, fac*x.^(1), 'r--', 'DisplayName', '$N$');
    end
end