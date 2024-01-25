# Data of OB22 (PPSN22)
Generate experimental results of OB22.
Run script (TRIALS can be changed accordingly), plots are generated automatically.
* fig1/main_fig1.m: generates left plot phi_i(sigma) and right sub-plots showing the error over all components
* fig2/main_fig2.m: R=100,10,1,0.1 is set at the top and generates sub-plots of Fig. 2.
* fig3/main_fig3_dyn.m: generate mean dynamics of ES. Choose population size of CONFIG 1 or CONFIG 2 (left/right sub-plots)
* fig3/main_fig3_iter.m: generate iterated dynamics using phi_i. Two configurations (1-2) by varying population size, and three sub-configurations (a-c) to generate dynamics using different phi-approximations (see Fig. 3)
(dynamics of real ES and iteration must be merged manually)
* fig4/main_fig4.m: 