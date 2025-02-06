# Data of: Bias in Standard Self-Adaptive Evolution Strategies.
Generate experimental results of OB24b (CEC 2024). 

Files to run simulation are contained in each fig-folder, respectively, and start with "main_":
* fig1/main_fig1.m: switch configurations A1, A2, B1, B2 to generate Rastrigin and Noisy Sphere dynamics using sigmaSA and CSA
* fig2/main_fig2.m: switch configurations 1 and 2 to generate dynamics of sigmaSA with log-normal and normal mutations
* fig3: main_contour.m creates progress rate landscape. main_dyn.m  is used to create corresponding dynamics of ES. The sigma*-R-dynamics is copied and pasted manually into progress rate landscape.
* fig4: main_vary.m: variable "SET_metaEP" sets (log-)normal mutations.
Variation is controlled by (un-)commenting data below DOE1 and DOE2, respectively.
A "*.mat" file is generated, which can then be evaluated using plots.m.
Data of log-normal and normal experiments are merged by copy and paste into a single plot.
* fig5/main_noise.m: config A creates dynamics using log-normal mutations, and config B the dynamics using normal mutations.
* fig6/main_cos.m: config A creates dynamics using log-normal mutations, and config B the dynamics using normal mutations.