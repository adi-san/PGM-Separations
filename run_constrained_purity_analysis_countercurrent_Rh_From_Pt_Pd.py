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



n_stages_int_list=[]
max_recov_list=[]
vol_flow_list=[]
C_out_list=[]
q_out_list=[]
max_sales_Rh_list=[]

price_arr=np.array([33730.32258,32929.67742, 164946.129])/1000*MW_arr # $/kg PGM --> $/mol PGM, Source, Johnson Matthey
# assume that these prices only apply in a particular stream rich with a component

Q_org_ig=.3

for n_stages_float in n_stages_arr:
  n_stages = int(n_stages_float)
  n_stages_int_list.append(n_stages)
  low_b = np.zeros(len(C_in)*n_stages)
  up_b=np.ones(len(C_in)*n_stages)*np.inf
  opt_bounds = (low_b, up_b)

  result=fsolve(Rh_purity_resid_fcn_countercurrent,Q_org_ig,args=(C_lig, Q_aq,n_stages, q_in,C_in,q_max_arr,K_Eq_arr,desired_purity_Rh), xtol=1e-6)
  
  Q_org_det=result[0]
#   print(result)
  print(f'For n_stages={n_stages}, the required organic flow rate to achieve {desired_purity_Rh}% purity of Rh in the aqueous phase is {Q_org_det} L/time.')
  C_counter_ig=np.tile(C_in,n_stages)
  C_countercurrent=least_squares(countercurrent_model, C_counter_ig, args=(C_lig, Q_aq, Q_org_det,n_stages, q_in,C_in,q_max_arr,K_Eq_arr), ftol=1e-11, gtol=1e-11, xtol=1e-11,bounds=opt_bounds, method='trf')
  C_countercurrent_concs=C_countercurrent.x
#   check the l2 norm of the objective function at the solution, which should be close to zero if the solution is good as mass is conserved across stages and thus the mass balance equations should be satisfied at the solution
  if np.linalg.norm(countercurrent_model(C_countercurrent_concs, C_lig, Q_aq, Q_org_det,n_stages, q_in,C_in,q_max_arr,K_Eq_arr))>1e-8:
    print(f'Warning: The solution for n_stages={n_stages} and Q_org={Q_org_det} may not have converged properly, as the l2 norm of the objective function is greater than 1e-8.')
  C_countercurrent_mat=C_countercurrent_concs.reshape(int(len(C_countercurrent_concs)/int(len(C_in))),int(len(C_in)))
  C_counter_plot=np.vstack((C_in,C_countercurrent_mat))
  recov_aq_arr=compute_recov_aq_arr(C_counter_plot,n_stages)
  Rh_recov=recov_aq_arr[-1,2]
  max_recov_list.append(Rh_recov)
  vol_flow_list.append(Q_org_det)

  # append the last stage exit conc vec to c_out_list as we will make a DF of these and plot them to see how the exit concentrations change with number of stages at the point where we achieve the desired purity value
  C_out_list.append(C_counter_plot[-1,:])
  q_out=q_func_langmuir_relation(q_max_arr, K_Eq_arr, C_counter_plot[1,:])
  q_out_list.append(q_out)
  price_Rh=price_arr[-1]*C_counter_plot[-1,-1] #$Rh/L feed
  max_sales_Rh_list.append(price_Rh)
  # print(price_Rh)

  # below is an approach to initialize guesses: simple bisection method
  if n_stages==2:
    Q_org_ig=Q_org_det/2
  elif n_stages==3:
    Q_org_ig=Q_org_det/1.5
  else:
    Q_org_ig=Q_org_det
  
print('Done with all simulations!')
max_recov_arr=np.array(max_recov_list)
vol_flow_arr=np.array(vol_flow_list)
#   solve the countercurrent model with the determined organic flow rate to get the recovery and purity values, which we can then plot to verify that we are indeed achieving the desired purity value and to see how the recovery changes with number of stages at the point of achieving the desired purity value. We can also plot the relationship between number of stages and required organic flow rate to achieve the desired purity value.
fig, ax1 = plt.subplots()
  # Recovery (left y-axis)
ax1.plot(n_stages_arr, max_recov_arr, linestyle='--', marker='v',color='tab:blue', label='Max Recovery')
ax1.set_xlabel('Stages in Countercurrent Operation')
ax1.set_ylabel('Max Recovery (%)', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')

  # Purity (right y-axis)
ax2 = ax1.twinx()
ax2.plot(n_stages_arr, vol_flow_arr, linestyle='--', marker='v',color='tab:red', label='Organic Flowrate')
# ax2.plot(test_flowrates, np.ones(len(test_flowrates))*95, color='tab:red',linestyle='--', label='Purity Threshold')
ax2.set_ylabel('Vol Flow Solution (L/time)', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')

  # Optional: combined legend
lines = ax1.get_lines() + ax2.get_lines()
labels: list[str] = [str(line.get_label()) for line in lines]
ax1.legend(lines, labels, loc='center right')
plt.title('Recovery at '+str(desired_purity_Rh)+'% Purity and Organic Flow Vs Stage Number')

plt.tight_layout()
plt.show()

# i want to append C_ to all the PGM_labels
Conc_out_labels=['C_'+label+'_out [M]' for label in PGM_labels]
Conc_in_labels=['C_'+label+'_in [M]' for label in PGM_labels]
q_in_labels=['q_'+label+'_in [mol PGM/mol ddFc]' for label in PGM_labels]
q_out_labels=['q_'+label+'_out [mol PGM/mol ligand]' for label in PGM_labels]
mat_bal_labels=['Material Residual '+label+' [mol PGM/time]' for label in PGM_labels]
# start generating a DF to save the results of the simulations in a nice format, which we can then save as a CSV and also use to make some nice plots to visualize how the exit concentrations, recoveries, and purities change with number of stages at the point where we achieve the desired purity value
stages_df=pd.DataFrame(n_stages_int_list, columns=['Stages'])
C_in_mat=np.tile(C_in,(highest_stage-starting_stage+1,1))
q_in_mat=np.tile(q_in,(highest_stage-starting_stage+1,1))
C_in_df=pd.DataFrame(C_in_mat, columns=Conc_in_labels)
q_in_df=pd.DataFrame(q_in_mat, columns=q_in_labels)

Q_org_arr_flat=np.tile(vol_flow_arr, len(C_in))
Q_org_mat=np.transpose(Q_org_arr_flat.reshape(len(C_in),int(len(Q_org_arr_flat)/len(C_in))))

C_out_mat=np.array(C_out_list)
q_out_mat=np.array(q_out_list)
C_out_df=pd.DataFrame(C_out_mat, columns=Conc_out_labels)
q_out_df=pd.DataFrame(q_out_mat, columns=q_out_labels)

mat_bal_mat=C_in_mat*Q_aq+q_in_mat*Q_org_mat*C_lig-C_out_mat*Q_aq-q_out_mat*Q_org_mat*C_lig
mat_bal_df=pd.DataFrame(mat_bal_mat, columns=mat_bal_labels)

result = pd.concat([stages_df,C_in_df, q_in_df, C_out_df, q_out_df, mat_bal_df], axis=1)
# result['Stages']=n_stages_int_list
result['Q_org [L/time]']=vol_flow_arr
result['Rh Recovery [%]']=max_recov_arr
result['PPM Mass Total [mg/L]']=total_conc_ppm_mass
result['Rh Relative Purity [%]']=desired_purity_Rh
result['Q_aq [L/time]']=Q_aq
result['C_lig [L/time]']=C_lig
result['Q_aq/(Q_org*C_lig) [1/M]']=Q_aq/(vol_flow_arr*C_lig)
result['Max Sales Rh [$Rh/L feed basis]']=max_sales_Rh_list
# result=result[['Stages']+Conc_in_labels+q_in_labels+Conc_out_labels+q_out_labels]
# result=result[['Stages']+Conc_in_labels+q_in_labels+Conc_out_labels+q
print(result)
result.to_csv('countercurrent_constrained_purity_analysis_results_'+str(total_conc_ppm_mass)+'_ppm_'+str(PGM_labels)+'_PGMs_'+str(desired_purity_Rh)+'_percent_purity_Rh_'+str(ligand)+'_ligand.csv', index=False)
# df.plot(linestyle='--',marker='v')

#   use the previous solution as the new guess, but scaled down by 10% to hopefully ensure that we approach the solution from the same direction each time and thus get a more monotonic relationship between number of stages and required organic flow rate
#   C_arr, q_arr = crosscurrent_model_fsolve(C_in, q_in, C_lig, Q_aq, Q_org_det, n_stages,q_max_arr,K_Eq_arr)
#   recov_aq_arr=compute_recov_aq_arr(C_arr,n_stages)
#   Rh_recov=recov_aq_arr[-1,2]
#   max_recov_list.append(Rh_recov)
#   vol_flow_list.append(Q_org_det)