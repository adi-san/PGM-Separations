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

q_max_label_list=['q_max_'+label for label in PGM_labels]
K_eq_label_list=['K_eq_'+label for label in PGM_labels]


# only fixed in first sep
rel_mol_frac_arr_S1=np.array([.45, .45, .1]) #Pt, Pd, Rh
C_in_S1=ppm_mass_tot_to_M_concs(total_conc_ppm_mass ,rel_mol_frac_arr_S1,MW_arr)
Q_aq_S1=1 # L/time
PGM_ligand_pairs_S1=[label+'_'+ligand_S1 for label in PGM_labels]
# print('PGM Ligand Pairs: ' +str(PGM_ligand_pairs))
PGM_ligand_pairs_S2=[label+'_'+ligand_S2 for label in PGM_labels]




Q_org_ig_S1=.2 #initial guess for organic flow rate, can be changed if needed

Q_org_ig_S2=.2 #initial guess for organic flow rate, can be changed if needed

starting_stage_S1=3
highest_stage_S1=12

starting_stage_S2=4
highest_stage_S2=13

q_max_arr_S1 = np.array([0.283184293825222,0.56945642916197,0.00376986289722275]) #mol PGM/mol ddFc
K_Eq_arr_S1 = np.array([1822.21079088769,2662.60590731049, 2401.50309384427]) # truly dimensionless parameters, 
# I may need to scale these by powers of 10 if convergence proves to be tricky

q_max_arr_S2 = np.array([0.018775015,0.104929611,0.001814625]) #mol PGM/mol coeFc
K_Eq_arr_S2 = np.array([285.1951444,2169.894489,866.4940906]) # truly dimensionless parameters, 
#I may need to scale these by powers of 10 if convergence proves to be tricky #coeFc

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
# print('Final Result DF from non-random functionalized code:')
# add K_eq and q_max values to the non_random_S1 df for clarity and so that I can use these values in the second separation
non_random_S1[q_max_label_list]=q_max_arr_S1
non_random_S1[K_eq_label_list]=K_Eq_arr_S1
# add S1 to all the columns of the non_random_S1 df to indicate that these results are from the first separation unit
non_random_S1.columns=['S1 '+col for col in non_random_S1.columns]
S1_column_list=non_random_S1.columns.to_list()
# print(non_random_S1)

non_random_df_list=[]
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
    non_random_S2[q_max_label_list]=q_max_arr_S2
    non_random_S2[K_eq_label_list]=K_Eq_arr_S2
    # print('Final Result DF from non-random functionalized code for second separation:')
    non_random_S2.columns=['S2 '+col for col in non_random_S2.columns]
    # print(non_random_S2)
    # create a combined df with the results from both separations
    # df1=non_random_S1.loc[index].to_frame().T.reset_index(drop=True)
    values_of_S1=row[S1_column_list].to_list()
    # make a df that is of the same dimensions as S2 but with the values of S1 repeated in each row
    df_S1_for_concat=pd.DataFrame([values_of_S1]*len(non_random_S2), columns=S1_column_list)
    combined_df=pd.concat([df_S1_for_concat.reset_index(drop=True), non_random_S2.reset_index(drop=True)], axis=1)
    # combined_df=pd.concat([non_random_S1.loc[index].to_frame().T.reset_index(drop=True), non_random_S2.reset_index(drop=True)], axis=1)
    non_random_df_list.append(combined_df)
    # print('Combined DF for index '+str(index)+':')
    # print(combined_df)
end=datetime.now()
print('Time taken for non-random case: '+str(end-start))
# concatenate the dfs
non_random_df=pd.concat(non_random_df_list,ignore_index=True)
# print(non_random_df)
non_random_df.to_csv('Constrained Purity Analysis Pt Pd Rh/constrained_purity_analysis_Pt_Pd_Rh_with_regressed_parameters.csv', index=False)

    # print('C_in_S2 for index '+str(index)+': '+str(C_in_S2))
    # get the C_in for the second separation from the results of the first separation
# we repeat the algorithm above for the randomly generated PGM parameters

df_isotherm_random=pd.read_csv('isotherm_parameters_with_random_samples.csv')

# I want to only look at the rows of this df that correspond to the same ligand and PGM pair (ex Pt_ddFc, Pd_ddFc, Rh_ddFc) as the non_random case, so I will filter the df based on the ligand and PGM_labels variables
# df_isotherm_random_filtered=df_isotherm_random[df_isotherm_random['PGM Ligand Pair']==PGM_ligand_pairs]
PGM_S1_ligand_pairs=[label+'_'+ligand_S1 for label in PGM_labels]

df_isotherm_random_filtered_S1 = df_isotherm_random[df_isotherm_random['PGM Ligand Pair'].isin(PGM_S1_ligand_pairs)]
# df_isotherm_random_filtered=df_isotherm_random_filtered[df_isotherm_random_filtered['PGM'].isin(PGM_labels)]
df_isotherm_random_filtered_q_max_S1=df_isotherm_random_filtered_S1[df_isotherm_random_filtered_S1['Parameter'].isin(['q_max [mol PGM/mol Ligand]'])]
df_isotherm_random_filtered_K_eq_S1=df_isotherm_random_filtered_S1[df_isotherm_random_filtered_S1['Parameter'].isin(['K_eq [~]'])]

# print(df_isotherm_random_filtered_q_max_S1)
PGM_S2_ligand_pairs=[label+'_'+ligand_S2 for label in PGM_labels]

df_isotherm_random_filtered_S2 = df_isotherm_random[df_isotherm_random['PGM Ligand Pair'].isin(PGM_S2_ligand_pairs)]
# df_isotherm_random_filtered=df_isotherm_random_filtered[df_isotherm_random_filtered['PGM'].isin(PGM_labels)]
df_isotherm_random_filtered_q_max_S2=df_isotherm_random_filtered_S2[df_isotherm_random_filtered_S2['Parameter'].isin(['q_max [mol PGM/mol Ligand]'])]
df_isotherm_random_filtered_K_eq_S2=df_isotherm_random_filtered_S2[df_isotherm_random_filtered_S2['Parameter'].isin(['K_eq [~]'])]
# for the whole shabang, do 500 random samples
num_samples=500
for i in range(num_samples):
    start=datetime.now()

    q_max_list_random_S1=df_isotherm_random_filtered_q_max_S1['Random Sample '+str(i+1)].to_list()
    K_Eq_list_random_S1=df_isotherm_random_filtered_K_eq_S1['Random Sample '+str(i+1)].to_list()
    q_max_arr_random_S1=np.array(q_max_list_random_S1)
    K_Eq_arr_random_S1=np.array(K_Eq_list_random_S1)

    q_max_list_random_S2=df_isotherm_random_filtered_q_max_S2['Random Sample '+str(i+1)].to_list()
    K_Eq_list_random_S2=df_isotherm_random_filtered_K_eq_S2['Random Sample '+str(i+1)].to_list()
    q_max_arr_random_S2=np.array(q_max_list_random_S2)
    K_Eq_arr_random_S2=np.array(K_Eq_list_random_S2)

    # print('q_max_arr_random for iteration '+str(i+1)+': '+str(q_max_arr_random))
    # print('K_Eq_arr_random for iteration '+str(i+1)+': '+str(K_Eq_arr_random))
    random_result_S1=run_constrained_purity_analysis_countercurrent(C_in_S1, 
                                                   q_in, 
                                                   C_lig, 
                                                   Q_aq_S1,
                                                   Q_org_ig_S1, 
                                                   q_max_arr_random_S1, 
                                                   K_Eq_arr_random_S1, 
                                                   desired_HK_PGM_purity,
                                                   starting_stage_S1,
                                                   highest_stage_S1,
                                                   ligand_S1,
                                                   PGM_labels,
                                                   MW_arr,
                                                   Rh_purity_resid_fcn_countercurrent)
    # print('Random Result DF from functionalized code for index '+str(index)+':')
    # q_max_K_eq_df=
    random_result_S1[q_max_label_list]=q_max_arr_random_S1
    random_result_S1[K_eq_label_list]=K_Eq_arr_random_S1

    random_result_S1.columns=['S1 '+col for col in random_result_S1.columns]
    S1_column_list=random_result_S1.columns.to_list()
    random_df_list=[]
    for index, row in random_result_S1.iterrows():
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
        random_S2=run_constrained_purity_analysis_countercurrent(C_in_S2, 
                                                    q_in, 
                                                    C_lig, 
                                                    Q_aq_S2,
                                                    Q_org_ig_S2, 
                                                    q_max_arr_random_S2, 
                                                    K_Eq_arr_random_S2, 
                                                    desired_HK_PGM_purity,
                                                    starting_stage_S2,
                                                    highest_stage_S2,
                                                    ligand_S2,
                                                    PGM_labels,
                                                    MW_arr,
                                                    Pt_purity_resid_fcn_countercurrent)
        random_S2[q_max_label_list]=q_max_arr_random_S2
        random_S2[K_eq_label_list]=K_Eq_arr_random_S2
        # print('Final Result DF from non-random functionalized code for second separation:')
        random_S2.columns=['S2 '+col for col in random_S2.columns]
        # print(non_random_S2)
        # create a combined df with the results from both separations
        # df1=non_random_S1.loc[index].to_frame().T.reset_index(drop=True)
        values_of_S1=row[S1_column_list].to_list()
        # make a df that is of the same dimensions as S2 but with the values of S1 repeated in each row
        df_S1_for_concat=pd.DataFrame([values_of_S1]*len(random_S2), columns=S1_column_list)
        combined_df=pd.concat([df_S1_for_concat.reset_index(drop=True), random_S2.reset_index(drop=True)], axis=1)
        # combined_df=pd.concat([non_random_S1.loc[index].to_frame().T.reset_index(drop=True), non_random_S2.reset_index(drop=True)], axis=1)
        random_df_list.append(combined_df)
    end=datetime.now()
    print('Time taken for iteration '+str(i+1)+': '+str(end-start))
    random_df=pd.concat(random_df_list,ignore_index=True)
    random_df.to_csv('Constrained Purity Analysis Pt Pd Rh/Uncertainty Results/constrained_purity_analysis_Pt_Pd_Rh_with_regressed_parameters_random_sample_'+str(i+1)+'.csv', index=False)
    end=datetime.now()
    # monte_carlo_df_list.append(random_result)

end_outer_loop=datetime.now()
