# Data of: Self-Adaptation of Multi-Recombinant Evolution Strategies on the Highly
Multimodal Rastrigin Function
Generate experimental results of OB24d (IEEE TEVC). 
Files to run simulation are contained in each fig-folder, respectively, and start with "main_":
* fig1/main_psi.m: run configs A1, A2, B1, B2 to generate respective sub-figures
* fig2/main_dyn.m: run configs A and B to generate dynamics. main_iter_sa.m: run configs A and B (swtich between 'R' and 'Y') to generate dynamics using progress rate and self-adaptation response, using both y- and R-iteration. Merge manually into dynamics of real ES.
* fig3/main_isolines_phi_psi.m: generate progress landscape for A=1 and A=7. Use main_dyn.m to generate median dynamics of real ES-runs. Use main_iter_sa.m to generate dynamics using R-iteration. Merge dynamics manually into progress landscape.
* fig4_1/main_isolines_phi_psi.m: generate progress landscape, main_iter_sa_terms.m: generate iteration (set NEGLECT_TERMS = 1 to neglect variance terms as shown in paper)
* fig4_2/main_isolines_phi_psi.m: generate progress landscape, use main_signcr.m to generate approximations of sigma*_crit. Copy and paste the data points back into progress landscape.
* fig5/main_isolines_phi_psi.m: generate progress landscape for the two configurations; run main_dyn.n for two configruations, and by varying TAU (comment and uncomment desired values)
* fig6/main_ert_varyAll_N50.m: run script to generate data for analysis of success rate P_S and expected runtime E_r. After data has been generated, run plots.m script.
* fig7/main_isolines_phi_psi.m: script generates progress rate and self-adaptation response landscape. Use main_dyn.m and main_iter_sa.m from fig3 to crease sigma*-R-dynamics.
* fig8/main_signcr.m: (un-)comment the configruations 1,2,3 to generate data. Change plot axis scale and labels to match the variable (MU, N, A) which was varied.
* fig9/visualize.m: run to generate plot.
* fig10/main_ert_varyAll_N200.m: run script to generate data for analysis of success rate P_S and expected runtime E_r (analogous to fig6). run plots.m to generate plots.