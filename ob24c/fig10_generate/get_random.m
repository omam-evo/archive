
function [fconf, init_stop, readme] = get_random(N)
    
        fconf = Fitness('Random',N,[]);

        init_stop.Y_0 = 1*ones(N,1);    % arbitrary (zero not shown in log-plot, use ones)
        init_stop.SIGMA_0 = 1;          % arbitrary
        init_stop.R_STOP = nan;
        init_stop.F_STOP = nan;
        init_stop.SIGMA_STOP = nan;
        init_stop.F_EVAL = nan; 
        init_stop.G = 1e3; % critical for RAM bc. of initialization in ES

        readme = {['N',num2str(N)]};

end