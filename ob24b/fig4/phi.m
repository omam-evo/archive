function [res] = phi(x, MU, N, OPTIONS, e_10, e_11, e_20)

    if any(strcmpi(OPTIONS,'full'))
        res = e_10 * x .* (1 + x.^2/(2*MU*N))./sqrt(1+x.^2/(2*N))./sqrt(1+x.^2/(MU*N)) - N*(sqrt(1+x.^2/(MU*N)) - 1);
    elseif any(strcmpi(OPTIONS,'medium'))
        res = e_10 * x ./ sqrt(1+x.^2/(2*N)) - x.^2/(2*MU);
    elseif any(strcmpi(OPTIONS,'large'))
        res = e_10*sqrt(2*N) - x.^2/(2*MU);
    elseif any(strcmpi(OPTIONS,'arn'))
        DQ = sqrt(1+x.^2/(2*N));
        res = x.*e_10./DQ - x.^2/(2*MU) .* ...
                (1 + 1/N*(e_11+(MU-1)*e_20)./(DQ.^2) - (N-1)/N^2.*x*e_10./DQ);
    elseif any(strcmpi(OPTIONS,'large_v2'))
        res = sqrt(2*N)*e_10 - e_20 - x.^2/(2*MU)*(1-sqrt(2/N)*e_10);
    else
        warning('Eq. type not defined.')
    end
    
end % fct