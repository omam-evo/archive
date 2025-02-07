# Data of: Mutation Strength Adaptation of the (µ/µI , λ)-ES for Large Population
Sizes on the Sphere Function

* fig1_left: main_csa_sa_ss.m generates dynamics on sphere for CSA and sigmaSA (see comments). Additionally, sigma^*_ss is measured and displayed together with N-dependent progress rate.
* fig1_right: main_progress_rates.m generates phi^* as a function of TAU. N-dependent phi^*, one-generation phi^* and large population approximation are shown.
* fig2: main_progress_rates.m creates Fig. 2a and 2b by (un-)commenting labeled sections. Increase TRIALS_PHI if needed.
* fig3: main_csa_iter.m generates iterated CSA dynamics of config "0" (MU=1000,N=100) in plot.
* fig4: main_csa_iter.m creates iterated steady-state results of Fig. 4a and 4b. Details given in header of file.
* fig5: see ob24c/fig2: the same steady-state CSA results are used in IEEE TEVC paper for population control strategies.
* fig6: main_sa.m runs steady-state simulations.
** LIST_ADAPT = {'N','2N','8N'} performs variation over TAU=1/sqrt(N),1/sqrt(2N),1/sqrt(8N) for each configuration.
** LIST_CONFIG = [1,2,3,4] performs variation over subfigures 1: (a), 2: (b), 3: (c), 4: (d)
** SET_META_EP = 0: lognormal mutations; SET_META_EP = 1: normal mutations (Meta-EP)
** FAC_TRIAL = 1: increase number of trials for smoother data.
** plots of log-normal and normal mutations need to be merged. Single plot contains data of all three TAU-values for each configuration.
* fig7: main_fig_sa.m: switch MODE = {'SA','METAEP'} to generate sphere dynamics with TAU = 1/sqrt(N). If necessary, change SEED=3 to different value to get different results.


