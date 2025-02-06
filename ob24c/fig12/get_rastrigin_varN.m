
function [fconf, init_stop, readme] = get_rastrigin_varN(MU, LAM)
    
    ALPHA = 2*pi;
    N_LIST = [10,30,100,300,1000];
    A = 3;

    for i=1:length(N_LIST)
        [fconf(i), init_stop(i), readme{i}] = get_rastrigin(N_LIST(i), A, ALPHA, MU, LAM);   
    end

end