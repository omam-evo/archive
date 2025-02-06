% Measuring # times of fitness deterioration
% L=10 means 10 differences evaluated, start evaluation at L+1
function [N_up, T] = pc_apop(pc_input, g, N_up, F_med_g)
    g0 = g-pc_input.L+1;  
    df = diff(F_med_g(g0:g));
    %% regular
    N_up(g+1) = sum(df>0)/(pc_input.L-1);
    %% weighted
    % w = 1:1:pc_input.L;
    % N_up(g+1) = sum(double(df>0).*w')/sum(w);

    x = N_up(g+1);
    if x < pc_input.apop_thresh
        T = 1;
    elseif x== pc_input.apop_thresh
        T = 0;
    else
        T = -1;
    end 
    % if g>=160
    %     dummy = 1;
    % end

end

