Source code for manuscript "Multiclass classification of growth curves using random change points and heterogeneous random effects"
by Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and Scott A. Sisson.

For questions, comments or remarks about the code please contact Vincent Chin.

/Code/

mcmc.m is the main file for the MCMC algorithm.

estimation.m is a sample script for running the MCMC algorithm on a simulated dataset (run data_simulation.m first in the /Simulation/Data/ folder).



/Simulation/Data/

data_simulation.m : MATLAB code for simulating data

label.mat : truth label (both df_fixed.mat and df_random.mat have the same truth labels)



/Simulation/Intermediate Results/

df_mf_*.Rdata : Results for fixed knot locations model fitted on data generated from fixed knot locations model

df_mr_*.Rdata : Results for random knot locations model fitted on data generated from fixed knot locations model

dr_mf_*.Rdata : Results for fixed knot locations model fitted on data generated from random knot locations model

dr_mr_*.Rdata : Results for random knot locations model fitted on data generated from random knot locations model

label.Rdata : Truth label (both df_fixed.mat and df_random.mat have the same truth labels)


To reproduce results presented in the manuscript, run plots.R file on the intermediate results.


