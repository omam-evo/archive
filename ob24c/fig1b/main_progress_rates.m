clear



MU_LIST = [10,30,100,300,1000,3000,10000];
N_LIST = 100*ones(1,length(MU_LIST));
RES_SIGN_OPT = nan*ones(1,length(MU_LIST)); %[6.30055	11.2006	18.55065 27.9007 42.330751 44.25085];
RES_SIGN_ZERO = nan*ones(1,length(MU_LIST)); %[12.16613	24.92011	47.9651	84.35109	154.70009	268.45009];     
TRIAL_LIST = [100,100,40,40];
MODE_MU_N = 1;
configstr = 'conf2';
TRIALS_PHI = 1e4;
THETA = 0.5;

tic
sph = Sphere;
f1 = figure; hold on;
for i=1:length(MU_LIST)

    MU = MU_LIST(i);
    LAM = MU/THETA;
    N = N_LIST(i);
    e_10 = e_vartheta_a_b(MU/LAM,1,0);
    e_11 = e_vartheta_a_b(MU/LAM,1,1);
    e_20 = e_vartheta_a_b(MU/LAM,2,0);

    [x_opt(i), x_2ndZero(i)] = sph.signorm_opt_sph(e_10, MU, LAM, N, {'full'});

    signzero_approx(i) = (8*N)^(1/4)*sqrt(e_10*MU);
    signzero_approx_v2(i) = sqrt(sqrt(2*N)*e_10*(2*MU)*(1-e_10/sqrt(2*N))/(1-sqrt(2/N)*e_10));
    
    x = linspace(0,1.1*x_2ndZero(i),1001);
    sign_opt = 6^(3/8)/3^(1/2)*(e_10*MU)^(1/4)*N^(3/8);
    
    phi_arn = phi(x, MU, N, {'arn'}, e_10, e_11, e_20);
    phi_full = phi(x, MU, N, {'full'}, e_10, e_11, e_20);
    phi_medium = phi(x, MU, N, {'medium'}, e_10, e_11, e_20);
    phi_large = phi(x, MU, N, {'large'}, e_10, e_11, e_20);
    phi_large_v2 = phi(x, MU, N, {'large_v2'}, e_10, e_11, e_20);
        
    SIGMA_LIST = linspace(0,1.1*x_2ndZero(i),10);
    SIGMA_LIST(1) = SIGMA_LIST(2)/2;
    [phi_R_mean(:,i), phi_R_stderr(:,i)] = fct_phi(MU, LAM, N, 0, 0, TRIALS_PHI, SIGMA_LIST);

        %plot(SIGMA_LIST, phi_mu1000_n100, 'ko')

        errorbar(SIGMA_LIST, phi_R_mean(:,i), phi_R_stderr(:,i), 'k.');
        plot(x, phi_full, 'k-')
        %plot(x, phi_arn, 'c-')
        %plot(x, phi_medium, 'r-.')
        plot(x, phi_large, 'k--')
        % plot(x, phi_large_v2, 'k--')
        % plot(signzero_approx, 0, 'rx');
        % xline(sign_opt, 'r--');
        % plot(signzero_approx_v2, 0, 'kx');
        % yline(0)

end
toc

f2 = figure; hold on;
plot(MU_LIST,x_2ndZero, 'r-')
plot(MU_LIST,x_opt, 'b-')
plot(MU_LIST,signzero_approx, 'k-.')
set(gca, 'XScale','Log')
set(gca, 'YScale','Log')
