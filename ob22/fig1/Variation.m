classdef Variation
    properties
    end
    
    methods   
        function[pos_list] = get_position_list(~, x0, x1, num)
            pos_list = linspace(x0, x1, num);
        end
        
        function[y_mat] = get_Y_ones(~, pos_list, N, ~)
            y_mat = ones(N,1).*pos_list;
        end 
        
        function[y_mat] = get_Y_zerosExcept1(~, pos_list, N, ~)
            y_mat = ones(N,1).*pos_list;
            y_mat(2:end,:) = 0;
        end 
        
        function[y_mat] = get_Y_from_R(~, R_LIST, N, SEED)
            rng(SEED);
            v = randn(N, length(R_LIST));
            v_scaled = repmat(R_LIST ./ vecnorm(v, 2, 1), N, 1);
            y_mat = v .* v_scaled;
        end      
    end % methods
end % class