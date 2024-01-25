
A = 10;
ALPHA = 2*pi;
N = 100;
VARTHETA = 0.5;
NUM_VAL = 1001;
R_end = 1*sqrt(N);
PLOT = 0;

MU_LIST = [10,50,100:100:1000];
%MU_LIST = [10,100,1000,10000];

RES_SIGN = nan*MU_LIST;
RES_R = nan*MU_LIST;
sigma0 = get_sigma_thresh(ALPHA, A);

warning('off')
for m=1:length(MU_LIST)
    
    MU = MU_LIST(m);
    fprintf('MU=%i\n',MU);
    LAM = MU/VARTHETA; 
    [s, r] = isolines_crit_point(MU, LAM, N, A, ALPHA, NUM_VAL, R_end, PLOT, 1, sigma0);
    RES_SIGN(m) = s;
    RES_R(m) = r;
end
warning('on')

figure; hold on;
plot(RES_SIGN, RES_R, 'k.');
get_ylim = ylim;
ylim(get_ylim);

sign=linspace(0,2*max(RES_SIGN),1001);
R_sigma0 = sigma0*N./sign;
plot(sign, R_sigma0, 'r--');
xlabel('$\sigma^*$');
ylabel('$R$')

myfigstyle(gcf,8,4,9,9);

xlim([0,45]);
xticks([0:10:50]);
ylim([0,7]);
yticks([0:2:7]);
% savefig('crit.fig');
% exportgraphics(gcf,'crit.pdf');