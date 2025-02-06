
function [fconf, init_stop, readme] = get_sphere(N, MU, LAM)
    
        fconf = Fitness('Sphere', N);

        init_stop.Y_0 = 1*ones(N,1); % convenience
        init_stop.R_STOP = nan;
        init_stop.F_STOP = 1e-12;    % fixed over all N  (same as Rastrigin) 
        init_stop.F_EVAL = nan;   % check budget
        init_stop.SIGMA_STOP = 1e-12;
        SIGMA_NORM = sqrt( sqrt(N*(8*e_vartheta_a_b(MU/LAM,1,0)^2*MU^2 + N)) - N);  % minimize initialization effects
        init_stop.SIGMA_0 = SIGMA_NORM*norm(init_stop.Y_0)/N;
        init_stop.G = 1e5; 
        readme = {['S',num2str(N)]};

end