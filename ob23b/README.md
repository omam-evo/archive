# Data of OB23b: Progress Analysis of a Multi-Recombinative Evolution Strategy on the
Highly Multimodal Rastrigin Function (Theoretical Computer Science)
Generate experimental results of OB23b.

Run script (TRIALS can be changed accordingly), plots are generated automatically.
* fig1/main_visualize.m: generate Rastrigin contour plot (2D, left) and regular plot (1D, right)
* fig2/main_dyn.m: evaluation of quality gain (two configurations N=10, N=100) and its normal approximation. Evaluation is done in muComLam_sSA.m and optimization is stopped after histogram was generated.
* fig3/main_phi.m: one-generation experiment phi(sigma). 6 Figures are generated in the end. Fig. 1-3 show left sub-plot (i=2) with approximations phi_i, phi_i(k_i), and phi_i(d_i) with different colors. Figs. 4-6 show right sub-plot (i=12) analogously.
* fig4_fig6/main_phi.m: one-generation experiment evaluates both first (Fig. 4) and second order progress rate (Fig. 6) for the same configuration, using two different R-initializations (left & right sub-plots, respectively). Hence, both figures (Fig. 4, Fig. 6) and sub-figures (left, right) are generated running the script once.
* fig5_fig7/main_phi.m: same procedure as for fig4_fig6, experiments using larger MU and N.
* fig8/main_dyn.m: generates and displays dynamics of 100 trials and corresponding median. Calls fct_iter.m in the end to evaluate iteration using first and second order progress rate. Two configurations (1-2) for left and right sub-plot.
* fig9/main_dyn.m: dynamics and iteration. variation of sigma* (choose values from comment), fct_iter.m is called automatically to evaluate iterated dynamics
* fig10/main_dyn.m: dynamics and iteration. The two configurations (1,2: left and right sub-plots) can be evaluated by (un)-commenting