
function [fconf, init_stop, readme] = get_random_var()
    
    N_LIST = [10,100,1000];
    for i=1:length(N_LIST)
        [fconf(i), init_stop(i), readme{i}] = get_random(N_LIST(i));
    end
end