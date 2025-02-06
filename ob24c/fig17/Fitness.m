classdef Fitness

    properties
        name
        N
        f
        y_hat
        label
        noise
        sigma_eps_g
        norm_noise
    end
    
    methods
        function obj = Fitness(name, N, params)
            obj.name = name;
            obj.N = N;
            if strcmpi(name, 'Sphere')
                obj.f = @(x) sum(x.^2, 1);
                obj.label = '$f_\mathrm{sph}$';
                obj.y_hat = zeros(N,1);
                obj.noise = 0;
            elseif strcmpi(name, 'Rastrigin')
                A = params(1);
                ALPHA = params(2);
                obj.label = '$f_\mathrm{ras}$';
                obj.f = @(x) sum(x.^2 + A*(1 - cos(ALPHA * x)), 1);
                obj.y_hat = zeros(N,1);
                obj.noise = 0;
            elseif strcmpi(name, 'SphereNoise')
                obj.f = @(x) sum(x.^2, 1);  % NOISE ADDED SEPARATELY
                obj.label = '$f_\mathrm{n}$';
                obj.y_hat = zeros(N,1);
                obj.noise = 1;
                obj.norm_noise = 0;
                obj.sigma_eps_g = params;
            elseif strcmpi(name, 'SphereNormNoise')
                obj.f = @(x) sum(x.^2, 1);  % NOISE ADDED SEPARATELY
                obj.label = '$f_\mathrm{nn}$';
                obj.y_hat = zeros(N,1);
                obj.noise = 1;
                obj.norm_noise = 1;
                obj.sigma_eps_g = params;
            elseif strcmpi(name, 'Linear')
                obj.label = '$f_\mathrm{lin}$';
                obj.f = @(x) sum(x,1);
                obj.y_hat = nan;
                obj.noise = 0;
            elseif strcmpi(name, 'Ellipsoid-1')
                obj.label = '$f_\mathrm{ell1}$';
                obj.f = @(x) (1:1:N) * (x.*x);    % dot-prod contains sum
                obj.y_hat = zeros(N,1);
                obj.noise = 0;
            elseif strcmpi(name, 'Ellipsoid-2')
                obj.label = '$f_\mathrm{ell2}$';
                obj.f = @(x) (1:1:N).^2 * (x.*x);    % dot-prod contains sum
                obj.y_hat = zeros(N,1);
                obj.noise = 0;
            elseif strcmpi(name, 'Ellipsoid-H')
                obj.label = '$f_\mathrm{ellH}$';
                obj.f = @(x) 10.^(6*((0:1:N-1))/(N-1)) * (x.*x);   % dot-prod contains sum
                obj.y_hat = zeros(N,1);
                obj.noise = 0;
            elseif strcmpi(name, 'EllipsoidHNoise')
                obj.label = '$f_\mathrm{nellH}$';
                obj.f = @(x) 10.^(6*((0:1:N-1))/(N-1)) * (x.*x);   % dot-prod contains sum
                obj.y_hat = zeros(N,1);
                obj.noise = 1;
                obj.norm_noise = 0;
                obj.sigma_eps_g = params;
            elseif strcmpi(name, 'Bohachevsky')
                obj.label = '$f_\mathrm{boh}$';
                obj.f = @(x) sum(params(1)+params(2)+x(1:size(x,1)-1,:).^2-params(1)*cos(params(3)*x(1:size(x,1)-1,:)) + ...
                    2*x(2:size(x,1),:).^2 - params(2)*cos(params(4)*x(2:size(x,1),:)) , 1);
%                 obj.f = @(y) 3 * sum(y.^2) ...
%                 - 2*y(1)^2 ...
%                 - y(N)^2 ...
%                 - 0.3 * (sum(cos(3 * pi * y))- cos(3 * pi * y(N))) ...
%                 - 0.4 * (sum(cos(4 * pi * y)) - cos(4 * pi * y(1))) ...
%                 + (N - 1)*(0.3 + 0.4);
                obj.y_hat = zeros(N,1);
                obj.noise = 0;
            elseif strcmpi(name, 'Random')
                obj.label = '$f_\mathrm{ran}$';
                obj.f = @(x) randn(1,size(x,2));
                obj.y_hat = zeros(N,1);
                obj.noise = 0;
            elseif strcmpi(name, 'Cosine')
                A = params(1);
                ALPHA = params(2);
                obj.label = '$f_\mathrm{cos}$';
                obj.f = @(x) sum(A*(1 - cos(ALPHA * x)), 1);
                obj.y_hat = zeros(N,1);
                obj.noise = 0;
            else
                error('Fitness not found')
            end
            
            
        end
%       function f = get_f(obj, x)
%           f = obj.f(x);
%       end
%       function f = get_fitness(obj, x, noise, R)
%             if nargin == 2 || noise==0  % prevent randn() from being called if no noise
%                 f = obj.f(x);
%             else
%                 if isnan(R) %R=nan => const noise strength
%                     f = obj.f(x) + noise*randn(1, size(x,2));
%                 else %R real valued => normalized noise strength
%                     f = obj.f(x) +(noise*2*R^2/size(x,1))*randn(1, size(x,2));
%                 end
%             end
%       end
      
    function f = get_f(obj, x, g, R)
        if obj.noise == 0 || nargin == 2
            f = obj.f(x);
        else
            if obj.norm_noise == 0
                f = obj.f(x) + obj.sigma_eps_g(g)*randn(1, size(x,2));
            else
                f = obj.f(x) + (obj.sigma_eps_g(g)*2*R^2/size(x,1))*randn(1, size(x,2));
            end
        end
      end     
      
      function [] = test(obj, x, y)
          if nargin == 2
            disp(num2str(x))
          else
            disp(num2str(x))
            disp(num2str(y))
          end
      end
        
    end
end

