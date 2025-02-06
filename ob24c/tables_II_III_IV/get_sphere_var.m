
function [fconf, init_stop, readme] = get_sphere_var(MU, LAM)

    N_LIST = [10,100,1000];
    for i=1:length(N_LIST)
        [fconf(i), init_stop(i), readme{i}] = get_sphere(N_LIST(i), MU, LAM);
    end

end