
MU_LIST = [10,100,1000,5000,1e4];
colors = {'r-','g-','b-','m-','k-'};
N = 100;
c = 0.9;
mode = 'medium';

figure; hold on;
for m=1:length(MU_LIST)
    MU = MU_LIST(m);
    LAM = MU/0.5;
    sph = Sphere;
    [x_opt, x_2ndZero] = sph.signorm_opt_sph(e_mu_lam_a_b_v2(MU,LAM,1,0), MU, LAM, N, {mode});

    x = linspace(0,x_2ndZero*1.1,1001);

    phi = sph.phi(x, e_mu_lam_a_b_v2(MU,LAM,1,0), MU, LAM, N, {mode});
    
    plot(x,phi, colors{m},'DisplayName',['mu=',num2str(MU)]);
    plot(x_2ndZero*c, sph.phi(x_2ndZero*c, e_mu_lam_a_b_v2(MU,LAM,1,0), MU, LAM, N, {mode}), 'ko');
    
    %% Predicted values
    signorm_pred = c*(8*N)^(1/4)*(e_vartheta_a_b(MU/LAM, 1, 0)*MU)^(1/2);
    phi_pred = sqrt(2*N)*e_vartheta_a_b(MU/LAM, 1, 0)*(1-c^2);
    plot(signorm_pred, phi_pred, 'rx');
    
end