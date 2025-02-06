clear
addpath('../n0_generic')

VARTHETA = 0.5;
e_10 = e_vartheta_a_b(VARTHETA, 1, 0);

%% SIGNLE
A_LIST = [1];
ALPHA_LIST = [2*pi];
N_LIST = [100];
MU_LIST = [100];
NUM_VAR = length(A_LIST);
VAR_LIST = MU_LIST;
SHOW_PLOT = 1;

%% MULTI

% A_LIST = [1,1,1,1,1,1,1]; %[0.5,1,2,4,6,8,10]; [1,1,1,1,1,1,1]
% ALPHA_LIST = [1*pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi];
% N_LIST = [100,100,100,100,100,100,100]; %[10,20,50,100,200,500,1000]; %[100,100,100,100,100,100,100];
% MU_LIST = [100,100,100,100,100,100,100]; %[10,30,50,100,300,500,1000]; [100,100,100,100,100,100,100];
% VAR_LIST = ALPHA_LIST;
% NUM_VAR = length(A_LIST);
% SHOW_PLOT = 0;

R_inf = nan*zeros(NUM_VAR,1); 
s0_iso = nan*zeros(NUM_VAR,1);
r0_iso = nan*zeros(NUM_VAR,1);
s0_phi = nan*zeros(NUM_VAR,1);
s0_w = nan*zeros(NUM_VAR,1);
s0_poly = nan*zeros(NUM_VAR,1);
s0_simpli = nan*zeros(NUM_VAR,1);
s_2ndZero = nan*zeros(NUM_VAR,1);

for i=1:NUM_VAR
    A = A_LIST(i);
    ALPHA = ALPHA_LIST(i);
    N = N_LIST(i);
    MU = MU_LIST(i);
    LAM = round(MU/VARTHETA);
    TAU = 1/sqrt(2*N);
    
    s_2ndZero(i) = sqrt( sqrt(N*(8*e_10^2*MU^2 + N)) - N);
    
    fprintf('A=%f, N=%i, MU=%i, ALPHA=%f\n',A,N,MU,ALPHA);
    [s_iso,r_iso] = func_isolines_phi_psi_newSS(MU,LAM,N,A,ALPHA,TAU,SHOW_PLOT);
    [sign0_phi, sign0_w, sign0_polyn1] = sign_crit(MU,LAM,N,A,ALPHA);
    R_inf(i) = sqrt(sqrt(N/2)*A*N/(4*e_10*MU));
    
    if SHOW_PLOT==1
        hold on;
        plot(s_iso,r_iso, 'r+')
        plot(sign0_phi,R_inf(i), 'bx')
        plot(sign0_w,R_inf(i), 'c*')
    end
   
    
    s0_iso(i) = s_iso;
    r0_iso(i) = r_iso;
    s0_phi(i) = sign0_phi;
    s0_w(i) = sign0_w;
    s0_poly(i) = sign0_polyn1;
    
    s0_simpli(i) = (512*N)^(1/4)*sqrt(e_10*MU)/ALPHA/sqrt(A)*sqrt(3)*sqrt(log(ALPHA*sqrt(A)/2));
    if imag(s0_simpli(i))~= 0
        s0_simpli(i) = nan;
    end
end

% figure; hold on; legend;
%     plot(VAR_LIST, s0_iso, 'ko-', 'DisplayName', 'Num.')
%     plot(VAR_LIST, s0_phi, 'co-','DisplayName', '$\varphi_R^{\mathrm{II},*}$=0')
%     plot(VAR_LIST, s0_w, 'mo-','DisplayName', '$W_0$')
%     plot(VAR_LIST, s0_poly, 'ro-','DisplayName', 'Polyn.')
%     plot(VAR_LIST, s_2ndZero, ':', 'color',[0.6,0.6,0.6], 'DisplayName', '$\sigma^*_{\varphi_0}$')
%     %plot(VAR_LIST, s0_simpli, 'o-','DisplayName', 'simpli')
% 
%     xlabel('$\mu N A$');
%     ylabel('$\sigma^*_{\mathrm{crit}}$');
% 
%     set(gca, 'XScale', 'log');
%     set(gca, 'YScale', 'log');
% 
%     myfigstyle(gcf, 16, 8, 10, 10);  

% figure; hold on; legend;
%     plot(VAR_LIST, r0_iso, 'ko-', 'DisplayName', 'R: contour');
%     plot(VAR_LIST, R_inf, 'o-','DisplayName', 'Rinf');

% factor_w = nan*A_LIS  T;
% for i=1:NUM_VAR
%     A = A_LIST(i);
%     factor_w(i) = lambertw(0, ALPHA^3*A^(3/2)/8);
% end
% figure;
%     plot(A_LIST, factor_w, 'k-')