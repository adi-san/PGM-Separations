from helper_functions import *

q_max_arr = np.array([0.283186498156831,0.569456410590181,0.00376984198947194]) #mol PGM/mol ddFc
K_Eq_arr = np.array([1822.14319447879,2662.60543887326, 2401.64491456341]) # truly dimensionless parameters, I may need to scale these by powers of 10 if convergence proves to be tricky

q_in=np.array([0,0,0])
total_conc=0.02 #units=M
C_in=np.array([.45, .45, .1])*total_conc
Q_aq=1 # L/time
C_lig=0.1 #mols of ligand/L solution

starting_stage=2
highest_stage=15
# the best number of stages is between 3 and 67
n_stages_arr=np.linspace(starting_stage,highest_stage,highest_stage-starting_stage+1)

# lower_bounds = np.array([0.05])
# upper_bounds = np.array([2.0])
# bounds = Bounds(lower_bounds, upper_bounds)

# specify this
desired_purity_Rh=95


max_recov_list=[]
vol_flow_list=[]

Q_org_ig=0.4



for n_stages_float in n_stages_arr:
  n_stages = int(n_stages_float)
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
  Q_org_ig=Q_org_det
print('Done with all simulations!')
max_recov_arr=np.array(max_recov_list)
vol_flow_arr=np.array(vol_flow_list)
#   solve the countercurrent model with the determined organic flow rate to get the recovery and purity values, which we can then plot to verify that we are indeed achieving the desired purity value and to see how the recovery changes with number of stages at the point of achieving the desired purity value. We can also plot the relationship between number of stages and required organic flow rate to achieve the desired purity value.

fig, ax1 = plt.subplots()

  # Recovery (left y-axis)
ax1.plot(n_stages_arr, max_recov_arr, color='tab:blue', label='Max Recovery')
ax1.set_xlabel('Stages in Countercurrent Operation')
ax1.set_ylabel('Max Recovery (%)', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')

  # Purity (right y-axis)
ax2 = ax1.twinx()
ax2.plot(n_stages_arr, vol_flow_arr, color='tab:red', label='Organic Flowrate')
# ax2.plot(test_flowrates, np.ones(len(test_flowrates))*95, color='tab:red',linestyle='--', label='Purity Threshold')
ax2.set_ylabel('Vol Flow Solution (L/time)', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')

  # Optional: combined legend
lines = ax1.get_lines() + ax2.get_lines()
labels: list[str] = [str(line.get_label()) for line in lines]
ax1.legend(lines, labels, loc='center right')
plt.title('Recovery at 95% Purity and Organic Flow Vs Stage Number')

plt.tight_layout()
plt.show()





  # use the previous solution as the new guess, but scaled down by 10% to hopefully ensure that we approach the solution from the same direction each time and thus get a more monotonic relationship between number of stages and required organic flow rate

#   C_arr, q_arr = crosscurrent_model_fsolve(C_in, q_in, C_lig, Q_aq, Q_org_det, n_stages,q_max_arr,K_Eq_arr)
#   recov_aq_arr=compute_recov_aq_arr(C_arr,n_stages)
#   Rh_recov=recov_aq_arr[-1,2]
#   max_recov_list.append(Rh_recov)
#   vol_flow_list.append(Q_org_det)