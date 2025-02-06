tau = @(A,ALPHA,N) sqrt(1/N.*(1-8./(ALPHA.^2.*A).*lambertw(0, ALPHA.^3.*A.^(3/2)/8)));

ALPHA = 2*pi;
A_LIST = linspace(0,10,1001);
N = 100;
res = tau(A_LIST,ALPHA,N);


% figure; hold on;
%     plot(A_LIST, res, 'k-')
%     yline(1/sqrt(N), 'r--')
%     yline(1/sqrt(2*N), 'b--')
%     xline(A_min)
    
x_min=1.47;
A_min = (2*x_min/ALPHA)^2;

x= 10.^linspace(-2,3,1001);
figure; hold on;
    plot(x, 2./x.^2.*lambertw(0, x.^3));
    %plot(x, 2./x.^2.*log((2*x.^3+1)./(1+log(x.^3+1))), 'r:');
    xline(x_min);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel('$x = \alpha\sqrt{A}/2$');
    ylabel('$2W_0(x^3)/x^2$');
    myfigstyle(gcf, 8, 4, 9, 9);
    xlim([1e-2,1e3]);
    xticks([1e-2,1e-1,1,1e1,1e2,1e3])
    yticks([1e-4,1e-2,1])