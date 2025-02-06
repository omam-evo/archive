
function [fconf, init_stop, readme] = get_rastrigin(N, A, ALPHA, MU, LAM)
    
        fconf = Fitness('Rastrigin', N, [A, ALPHA]);
    
        init_stop.Y_0 = 2*ceil(ALPHA*A/2)*ones(N,1); %slightly larger than last attractor
        init_stop.R_STOP = nan;
        init_stop.F_STOP = 1e-3;  % in GA optimizing sphere (to some extent)
        init_stop.F_EVAL = 5e7;   % limit of resources
        init_stop.SIGMA_STOP = 1e-5;
        SIGMA_NORM = sqrt( sqrt(N*(8*e_vartheta_a_b(MU/LAM,1,0)^2*MU^2 + N)) - N); % minimize initialization effects
        init_stop.SIGMA_0 = SIGMA_NORM*norm(init_stop.Y_0)/N;
        init_stop.G = 1e5; 
        readme = {['Ras.~$N$=',num2str(N), ' $A$=',num2str(A), ' $\alpha$=',num2str(ALPHA/pi),'$\pi$']};
end