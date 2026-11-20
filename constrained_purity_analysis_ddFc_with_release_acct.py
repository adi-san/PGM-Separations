from helper_functions import *

# this is an example of constrained analysis for a single ligand

ligand='ddFc'
PGM_labels=['Pt','Pd','Rh']
q_max_arr = np.array([0.283184293825222,0.56945642916197,0.00376986289722275]) #mol PGM/mol ddFc
# q_max_arr = np.array([0.018775015,0.104929611,0.001814625]) #mol PGM/mol coeFc
# q_max_arr = np.array([0.272571301,0.409453841,0.003109267]) #mol PGM/mol deFc

K_Eq_arr = np.array([1822.21079088769,2662.60590731049, 2401.50309384427]) # truly dimensionless parameters, I may need to scale these by powers of 10 if convergence proves to be tricky #ddFc
# K_Eq_arr = np.array([285.1951444,2169.894489,866.4940906]) # truly dimensionless parameters, I may need to scale these by powers of 10 if convergence proves to be tricky #coeFc
# K_Eq_arr = np.array([349.540848,950.3304042, 368.9954269]) # truly dimensionless parameters, I may need to scale these by powers of 10 if convergence proves to be tricky #deFc

# will need to make this something we solve for
q_in=np.array([0,0,0])

total_conc_ppm_mass=500 # mg/L
rel_mol_frac_arr=np.array([.45, .45, .1]) #Pt, Pd, Rh
MW_arr=np.array([195.084,106.42,102.91])
C_in=ppm_mass_tot_to_M_concs(total_conc_ppm_mass ,rel_mol_frac_arr,MW_arr)


C_lig=0.1 #mols of ligand/L solution
Q_aq=1 # L/time

starting_stage=3
highest_stage=12
# the best number of stages is between 3 and 67
n_stages_arr=np.linspace(starting_stage,highest_stage,highest_stage-starting_stage+1)

# specify this
desired_purity_Rh=95 # in percent, so 95 means 95% purity of Rh in the aqueous phase

# main change here is that we account for release effienccy of ligand and purge fraction
release_efficiency=0.4
purge_fraction=0.02


n_stages_int_list=[]
max_recov_list=[]
vol_flow_list=[]
C_out_list=[]
q_out_list=[]
max_sales_Rh_list=[]

price_arr=np.array([33730.32258,32929.67742, 164946.129])/1000*MW_arr # $/kg PGM --> $/mol PGM, Source, Johnson Matthey
# assume that these prices only apply in a particular stream rich with a component

Q_org_ig=.3

# for n_stages_float in n_stages_arr:
#     n_stages_int=int(n_stages_float)
#     n_stages_int_list.append(n_stages_int)
#     C_out,q_out,max_recov,vol_flow,max_sales_Rh=constrained_purity_analysis_ddFc_with_release_acct(C_in,C_lig,Q_aq,n_stages_int,PGM_labels,q_max_arr,K_Eq_arr,desired_purity_Rh,release_efficiency,purge_fraction,price_arr,Q_org_ig)
#     max_recov_list.append(max_recov)
#     vol_flow_list.append(vol_flow)
#     C_out_list.append(C_out)
#     q_out_list.append(q_out)
#     max_sales_Rh_list.append(max_sales_Rh)