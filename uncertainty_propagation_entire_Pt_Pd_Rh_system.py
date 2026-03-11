from helper_functions import *

# so we are going to propagate uncertainty throughout the first and second separation unit
# first unit will be ddFc ligand to separate out Rh from Pt and Pd
# second unit will be coeFc ligand to separate out Pt from Pd and Rh
# We will vary first unit stages from 3-12
# We will vary second unit stages from 4-13
# therefore, there are 10x10=100 possible pairs of stages
'''
* we will do 500 random generations of 12 parameters (q_i_maxs and K_eqs):
* for each of the 500 randomly generated combos:
    * run the constrained purity analysis for the first separation
    * these results will be in the form of a df
    * use this df and iterate thru each of the results, running the purity analysis for each of the first DF cases
'''

ligand_S1='ddFc'
ligand_S2='coeFc'
PGM_labels=['Pt','Pd','Rh']



# fixed in both sep units
q_in=np.array([0,0,0])
total_conc_ppm_mass=500 # mg/L
MW_arr=np.array([195.084,106.42,102.91])
C_lig=0.1 #mols of ligand/L solution
desired_HK_PGM_purity=95

# only fixed in first sep
rel_mol_frac_arr_S1=np.array([.45, .45, .1]) #Pt, Pd, Rh
C_in_S1=ppm_mass_tot_to_M_concs(total_conc_ppm_mass ,rel_mol_frac_arr_S1,MW_arr)
Q_aq_S1=1 # L/time
PGM_ligand_pairs_S1=[label+'_'+ligand_S1 for label in PGM_labels]
# print('PGM Ligand Pairs: ' +str(PGM_ligand_pairs))
PGM_ligand_pairs_S2=[label+'_'+ligand_S2 for label in PGM_labels]




Q_org_ig_S1=.2 #initial guess for organic flow rate, can be changed if needed

Q_org_ig_S2=.15 #initial guess for organic flow rate, can be changed if needed

starting_stage_S1=3
highest_stage_S1=12

starting_stage_S2=4
highest_stage_S2=13

q_max_arr_S1 = np.array([0.283184293825222,0.56945642916197,0.00376986289722275]) #mol PGM/mol ddFc
K_Eq_arr_S1 = np.array([1822.21079088769,2662.60590731049, 2401.50309384427]) # truly dimensionless parameters, I may need to scale these by powers of 10 if convergence proves to be tricky

q_max_arr_S2 = np.array([0.018775015,0.104929611,0.001814625]) #mol PGM/mol coeFc
K_Eq_arr_S2 = np.array([285.1951444,2169.894489,866.4940906]) # truly dimensionless parameters, I may need to scale these by powers of 10 if convergence proves to be tricky #coeFc

start=datetime.now()
non_random_S1=run_constrained_purity_analysis_countercurrent(C_in_S1, 
                                                   q_in, 
                                                   C_lig, 
                                                   Q_aq_S1,
                                                   Q_org_ig_S1, 
                                                   q_max_arr_S1, 
                                                   K_Eq_arr_S1, 
                                                   desired_HK_PGM_purity,
                                                   starting_stage_S1,
                                                   highest_stage_S1,
                                                   ligand_S1,
                                                   PGM_labels,
                                                   MW_arr,
                                                   Rh_purity_resid_fcn_countercurrent)
print('Final Result DF from non-random functionalized code:')
# add S1 to all the columns of the non_random_S1 df to indicate that these results are from the first separation unit
non_random_S1.columns=['S1 '+col for col in non_random_S1.columns]
print(non_random_S1)

for index, row in non_random_S1.iterrows():
    # run the second separation for each row of the first separation results df
    # print('Running second separation for index '+str(index)+' of the first separation results df...')
    # get the C_in for the second separation from the results of the first separation
    q_out_S1=row[['S1 q_'+label+'_out [mol PGM/mol ligand]' for label in PGM_labels]].to_numpy()
    Q_org_out_S1=row['S1 Q_org [L/time]']
    rel_mol_frac_arr_S2=q_out_S1/np.sum(q_out_S1) #Pt, Pd, Rh
    C_in_S2=ppm_mass_tot_to_M_concs(total_conc_ppm_mass ,rel_mol_frac_arr_S2,MW_arr)
    Q_aq_S2=Q_aq=water_flow_successive_stages(q_out_S1, 
                                 total_conc_ppm_mass, 
                                 MW_arr, 
                                 Q_org_out_S1,
                                 C_lig) # L/time
    non_random_S2=run_constrained_purity_analysis_countercurrent(C_in_S2, 
                                                   q_in, 
                                                   C_lig, 
                                                   Q_aq_S2,
                                                   Q_org_ig_S2, 
                                                   q_max_arr_S2, 
                                                   K_Eq_arr_S2, 
                                                   desired_HK_PGM_purity,
                                                   starting_stage_S2,
                                                   highest_stage_S2,
                                                   ligand_S2,
                                                   PGM_labels,
                                                   MW_arr,
                                                   Pt_purity_resid_fcn_countercurrent)
    # print('Final Result DF from non-random functionalized code for second separation:')
    non_random_S2.columns=['S2 '+col for col in non_random_S2.columns]
    print(non_random_S2)
end=datetime.now()
print('Time taken for non-random case: '+str(end-start))
print('mom')

    # print('C_in_S2 for index '+str(index)+': '+str(C_in_S2))
    # get the C_in for the second separation from the results of the first separation
