helper_functions.py contains the boilerplate functions used for nonlinear regression, uncertaint analysis, and countercurrent eLLE simulations.
run extract_isotherm_parameters.py to extract q_max,i and K_Eq,i parameters from experimental data along with their respective standard errors, approximated using a first order Taylor Series Approxiamtion.
generate_random_parameters can be used to create 500 random combinations of q_max,i and K_Eq,i parameters used for the Monte Carlo uncertainty analysis.
the "run_constrained_purity_analysis....py" files are used to generate the LLE balances for a given ligand, feed composition, and stage number with the Pd from Pt and Rh being the second separation and the Rh from Pt and Pd being the first separation.
The uncertainty_propagation_entire_Pt_Pd_Rh_system.py file is used to generate the LLE balances for a Pt-Pd-Rh system assuming the aqueous feed to each LLE column is 500 ppm PGM and also generates the Monte Carlo uncertainty analysis LLE material balances.
  THIS FILE WILL TAKE A FEW HOURS TO RUN, so the resultrs can be accessed under the Constrained Purity Analysis Pt Pd Rh folder.
Pt_Pd_Rh_S1_S2_data_cleaning_and_plotting.py is used to clean the Mante Carlo data, removing random parameter cases that generate physically unmeaningful results.
plot_uncertainty_Pt_Pd_Rh_system_S1_S2.py is used to generate the uncertainty figure where S1 stages vary on the x-axis.
run_countercurrent_simulations.py contains the code used to generate the countercurrent stage concentration profile plot.

