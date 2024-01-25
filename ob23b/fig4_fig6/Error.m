classdef Error 
    methods     (Static=true)
        
        %%
        % Error: sorted values (w.r.t. id_sorted)
        %        euclidian norm over all Sigma-Values for Position-Vector p and component n (n,p)
        function [error, id_sorted, phi_squared_error] = mean_squared(phi_meas, phi_approx)

            phi_squared_error = (phi_approx-phi_meas).^2;

            mean_squared = mean(phi_squared_error, 2); % Larger mean squared error if overall phi-niveau is larger
            relative_squared = sum(phi_squared_error, 2)./ sum((mean(phi_meas, 2)-phi_meas).^2, 2); % normalized by mean(phi_mean) values (phi-niveau)

            [N, NUM_SIG, NUM_POS] = size(phi_squared_error);
            %     if dim==1
            %         delta_phi = permute(delta_phi, [2,3,1]);
            %     elseif dim==2
            %         delta_phi = permute(delta_phi, [1,3,2]);
            %     end

            %% Index to be iterated over is last
            error = nan*zeros(N, NUM_POS);
            for p=1:NUM_POS
                error(:, p) = mean_squared(:,1,p); %vecnorm(delta_phi(:,:,p), 2, 2);
            end
            [~, id_sorted] = sort(error, 1, 'ascend');

        end
        function [error] = difference_at_location(phi_meas, phi_approx)
            error = phi_approx - phi_meas; 
        end
        
        function [error] = difference_at_location_relative(phi_meas, phi_approx, idx)
            error = (phi_approx(:, idx, :) - phi_meas(:, idx, :))./phi_meas(:, idx, :); 
        end
        
        %% GET error at max value and close to zero
        function [error_max,error_secZero] = difference_at_maxSigma_secZero(phi_meas, phi_approx)
            id = 1; 
            NUM_POS = size(phi_meas,2);
            [phi_max, id_max] = max(phi_meas, [], id);  % get max. phi over sigma
            error_max = nan*zeros(1,NUM_POS);           
            error_secZero = nan*zeros(1,NUM_POS);
            
            % Evaluate for each position p (all sigma values)
            for p=1:NUM_POS
                error_max(1,p) = phi_approx(id_max(p),p) - phi_max(1,p); % get error at max. location
                phi_after_max = phi_meas(id_max(p):end,p);               % get all phi values after maximum (2nd zero)
                [~, minId] = min(abs(phi_after_max));                    % find id of value closest to 0 (rel. to maximum location)
                minVal = phi_after_max(minId);                           % get actual value => can be larger/smaller 0
                error_secZero(1,p) = phi_approx(id_max(p)+minId-1,p) - minVal;  % get error of value clostest to 2nd zero
            end 
        end
        
        %% GET error at location given by id
        function [error_id] = difference_at_constSigma(phi_meas, phi_approx, id)
            error_id = phi_approx(id,:) - phi_meas(id,:);
        end
  
    end
end


