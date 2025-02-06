
function [fconf, init_stop, readme] = get_rastrigin_varAlpha(MU, LAM, omitconfig)
    
    A = 3;
    N_LIST = [10,30,100,300,1000];
    ALPHA_LIST = [10*pi,7*pi,4*pi,3*pi,2*pi];

    ids = 1:1:length(N_LIST);
    N_LIST = N_LIST(ids~=omitconfig);
    ALPHA_LIST = ALPHA_LIST(ids~=omitconfig);

    for i=1:length(N_LIST)
        [fconf(i), init_stop(i), readme{i}] = get_rastrigin(N_LIST(i), A, ALPHA_LIST(i), MU, LAM);
        
    end

end