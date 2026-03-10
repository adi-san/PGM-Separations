from helper_functions import *



ligand='ddFc'
PGM_labels=['Pt','Pd','Rh']

PGM_ligand_pairs=[label+'_'+ligand for label in PGM_labels]
# print('PGM Ligand Pairs: ' +str(PGM_ligand_pairs))

q_in=np.array([0,0,0])

total_conc_ppm_mass=500 # mg/L
rel_mol_frac_arr=np.array([.45, .45, .1]) #Pt, Pd, Rh
MW_arr=np.array([195.084,106.42,102.91])
C_in=ppm_mass_tot_to_M_concs(total_conc_ppm_mass ,rel_mol_frac_arr,MW_arr)

C_lig=0.1 #mols of ligand/L solution
Q_aq=1 # L/time

Q_org_ig=.2 #initial guess for organic flow rate, can be changed if needed
starting_stage=3
highest_stage=12

q_max_arr = np.array([0.283184293825222,0.56945642916197,0.00376986289722275]) #mol PGM/mol ddFc
K_Eq_arr = np.array([1822.21079088769,2662.60590731049, 2401.50309384427]) # truly dimensionless parameters, I may need to scale these by powers of 10 if convergence proves to be tricky #ddFc


df_isotherm_random=pd.read_csv('isotherm_parameters_with_random_samples.csv')
# I want to only look at the rows of this df that correspond to the same ligand and PGM pair (ex Pt_ddFc, Pd_ddFc, Rh_ddFc) as the non_random case, so I will filter the df based on the ligand and PGM_labels variables
# df_isotherm_random_filtered=df_isotherm_random[df_isotherm_random['PGM Ligand Pair']==PGM_ligand_pairs]

df_isotherm_random_filtered = df_isotherm_random[df_isotherm_random['PGM Ligand Pair'].isin(PGM_ligand_pairs)]
# df_isotherm_random_filtered=df_isotherm_random_filtered[df_isotherm_random_filtered['PGM'].isin(PGM_labels)]
df_isotherm_random_filtered_q_max=df_isotherm_random_filtered[df_isotherm_random_filtered['Parameter'].isin(['q_max [mol PGM/mol Ligand]'])]
df_isotherm_random_filtered_K_eq=df_isotherm_random_filtered[df_isotherm_random_filtered['Parameter'].isin(['K_eq [~]'])]

q_max_label_list=['q_max_'+label for label in PGM_labels]
K_eq_label_list=['K_eq_'+label for label in PGM_labels]

print(df_isotherm_random_filtered_q_max)
start_outer_loop=datetime.now()
monte_carlo_df_list=[]

num_samples=100
# takes about 4 hrs to run this with 1000 samples
for i in range(num_samples):
    start=datetime.now()
    q_max_list_random=df_isotherm_random_filtered_q_max['Random Sample '+str(i+1)].to_list()
    K_Eq_list_random=df_isotherm_random_filtered_K_eq['Random Sample '+str(i+1)].to_list()
    q_max_arr_random=np.array(q_max_list_random)
    K_Eq_arr_random=np.array(K_Eq_list_random)
    print('q_max_arr_random for iteration '+str(i+1)+': '+str(q_max_arr_random))
    print('K_Eq_arr_random for iteration '+str(i+1)+': '+str(K_Eq_arr_random))
    random_result=run_constrained_purity_analysis_countercurrent(C_in, 
                                                   q_in, 
                                                   C_lig, 
                                                   Q_aq,
                                                   Q_org_ig, 
                                                   q_max_arr_random, 
                                                   K_Eq_arr_random, 
                                                   95,
                                                   starting_stage,
                                                   highest_stage,
                                                   ligand,
                                                   PGM_labels,
                                                   MW_arr,
                                                   Rh_purity_resid_fcn_countercurrent)
    # print('Random Result DF from functionalized code for index '+str(index)+':')
    # q_max_K_eq_df=
    random_result[q_max_label_list]=q_max_arr_random
    random_result[K_eq_label_list]=K_Eq_arr_random
    print(random_result)
    end=datetime.now()
    print('Time taken for iteration '+str(i+1)+': '+str(end-start))
    monte_carlo_df_list.append(random_result)

end_outer_loop=datetime.now()
print('Total time taken for all iterations: '+str(end_outer_loop-start_outer_loop))

monte_carlo_df=pd.concat(monte_carlo_df_list, ignore_index=True)
print('Final Monte Carlo DF:')
print(monte_carlo_df)
monte_carlo_df.to_csv('monte_carlo_results_'+ligand+'_'+str(num_samples)+'_random_samples'+'_Rh_From_Pt_Pd'+'.csv', index=False)

non_random=run_constrained_purity_analysis_countercurrent(C_in, 
                                                   q_in, 
                                                   C_lig, 
                                                   Q_aq,
                                                   Q_org_ig, 
                                                   q_max_arr, 
                                                   K_Eq_arr, 
                                                   95,
                                                   starting_stage,
                                                   highest_stage,
                                                   ligand,
                                                   PGM_labels,
                                                   MW_arr,
                                                   Rh_purity_resid_fcn_countercurrent)
print('Final Result DF from non-random functionalized code:')
print(non_random)

