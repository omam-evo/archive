
function [fconf, init_stop, readme] = get_rastrigin_varA(MU, LAM, omitconfig)
    
    ALPHA = 2*pi;
    N_LIST = [10,30,100,300,1000];
    A_LIST = [65,33,12,7,3];

    ids = 1:1:length(N_LIST);
    N_LIST = N_LIST(ids~=omitconfig);
    A_LIST = A_LIST(ids~=omitconfig);

    for i=1:length(N_LIST)
        [fconf(i), init_stop(i), readme{i}] = get_rastrigin(N_LIST(i), A_LIST(i), ALPHA, MU, LAM);
        
    end

end