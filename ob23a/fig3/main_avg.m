%% Rastrigin as a function of R: f(y)->f(R)

N_LIST = [3,10,100,1000];
A = 10; ALPHA = 2*pi;
TR = 10^4;
ras = RasTools;

R_LIST = 10.^linspace(-1,3,101)';
NUM_VAR = length(R_LIST);
F_v1 = nan*zeros(NUM_VAR, 1);
F_v2 = nan*zeros(NUM_VAR, 1);
F_v3 = nan*zeros(NUM_VAR, 1);
F_v4 = nan*zeros(NUM_VAR, 1);
F_v5 = nan*zeros(NUM_VAR, 1);

%% Exact analytic averages
F_N1 = @(ALPHA,A,N,R) R^2 + N*A*(1 - cos(ALPHA*R));
F_bessel_N2 = @(ALPHA,A,N,R) R.^2 + 2*A * (1 - besselj(0,ALPHA*R));

%% Expected value average
F_EV = @(ALPHA,A,N,R) R^2 + N*A*(1 - exp(-0.5*(ALPHA*R)^2/N));

%%
rng(1);
fig=figure;hold on;
xlabel('$R$'); ylabel('$f(R)$'); 
%title(['$\alpha=2\pi$, $A=10$; bottom to top: $N=\{3,10,100,1000\}$']); 

for N = N_LIST
    fprintf('\t N=%i\n',N)
    for i=1:NUM_VAR
        R = R_LIST(i);

        %% v1: norm(Y)=R
        Y = randn(N, TR);
        Y = Y./vecnorm(Y, 2, 1)*R;
        F_v1(i) = mean(ras.f(Y, ALPHA, A));

        %% v2: Y~R/sqrt(N)*N(0,1)
        %Y = R/sqrt(N)*randn(N, TR);
        %F_v2(i) = mean(ras.f(Y, ALPHA, A));

        %% v3: R dependent in expectation
        if N==1
            F_v3(i) = R^2 + N*A*(1 - cos(ALPHA*R));  %R^2 + N*A*(1 - exp(-(ALPHA*R/sqrt(N))^2/2)); 
        end
        
        %% v4:
        if N==2
            F_v4(i) = F_bessel_N2(ALPHA,A,N,R);
        end
        
        if N>=3
            F_v5(i) = F_EV(ALPHA,A,N,R);
        end
               
    end
    plot(R_LIST, F_v1, 'k-', 'DisplayName', 'Sampling $||\textbf{\emph y}||=R$');
    %plot(R_LIST, F_v2, '--','DisplayName', '$y_i \sim R/\sqrt{N}\mathcal{N}(0,1)$');
    if N==1, plot(R_LIST, F_v3, 'g:','DisplayName', '$F(R)=R^2+A(1-\cos(\alpha R))$'); end
    if N==2, plot(R_LIST, F_v4, 'g:','DisplayName', '$F(R)=R^2+2A(1-J_0(\alpha R))$'); end
    if N>=3, plot(R_LIST, F_v5, 'g:','DisplayName', '$F(R)=R^2+NA(1-\mathrm{exp}[-(\alpha R)^2/N])$'); end
end

% C = 1; % transition constant
% for N = N_LIST
%     R_trans = sqrt(C*2*N)/ALPHA;
%     plot(R_trans, F_EV(ALPHA,A,N,R_trans), 'rx');
%     fprintf('N=%i, R_trans=%d \n',N,R_trans);
% end

%legend('Location', 'best');
xlim([R_LIST(1),R_LIST(end)]);
ylim([1e0,1e6]);
set(findall(gcf, '-property', 'XScale'), 'XScale', 'log');
set(findall(gcf, '-property', 'YScale'), 'YScale', 'log');
myfigstyle(gcf, 7.5, 4.5, 9, 9);
xticks([1e-1,1,1e1,1e2,1e3])
