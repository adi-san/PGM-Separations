from helper_functions import *

df_constrained_nonrandom=pd.read_csv('countercurrent_constrained_purity_analysis_results_500_ppm_[\'Pt\', \'Pd\', \'Rh\']_PGMs_95_percent_purity_Rh_ddFc_ligand.csv')

# print(df_constrained_nonrandom)

df_random=pd.read_csv('monte_carlo_results_ddFc_Rh_From_Pt_Pd.csv')
# df_random=pd.read_csv('monte_carlo_results_ddFc_100_random_samples_Rh_From_Pt_Pd.csv')
print(df_random)
# in the random df, I want to remove any of the rows where q_org is negative or greater than 1, as these indicate nonsense results
df_random_filtered=df_random[(df_random['Q_org [L/time]']>0) & (df_random['Q_org [L/time]']<1)]
# get the K_eq and q_max values from the rows we dropped, and see if they are in a reasonable range
dropped_rows=df_random[~((df_random['Q_org [L/time]']>0) & (df_random['Q_org [L/time]']<1))]
print('Dropped rows with invalid Q_org values:')
print(dropped_rows[['K_eq_Pt', 'q_max_Pt', 'K_eq_Pd', 'q_max_Pd', 'K_eq_Rh', 'q_max_Rh']])
# in df_random_filtered, i want to drop the rows containing any of the combos of K_eq and q_max values that are in the dropped_rows, since those combos seem to lead to invalid results
for index, row in dropped_rows.iterrows():
    K_eq_Pt=row['K_eq_Pt']
    q_max_Pt=row['q_max_Pt']
    K_eq_Pd=row['K_eq_Pd']
    q_max_Pd=row['q_max_Pd']
    K_eq_Rh=row['K_eq_Rh']
    q_max_Rh=row['q_max_Rh']
    df_random_filtered=df_random_filtered[~((df_random_filtered['K_eq_Pt']==K_eq_Pt) & (df_random_filtered['q_max_Pt']==q_max_Pt) & (df_random_filtered['K_eq_Pd']==K_eq_Pd) & (df_random_filtered['q_max_Pd']==q_max_Pd) & (df_random_filtered['K_eq_Rh']==K_eq_Rh) & (df_random_filtered['q_max_Rh']==q_max_Rh))]

# print(df_random_filtered)
# sort the df_random_filtered by the number of stages
# df_random_filtered_sorted=df_random_filtered.sort_values(by='Stages')
# for each chunk where each chunk corresponds to a different number of stages, I want the min and max recovery nad the min and max organic flow rate

# count=0
stage_list=[]
min_recovery_list=[]
max_recovery_list=[]
min_flow_list=[]
max_flow_list=[]
for stages, group in df_random_filtered.groupby('Stages'):
    min_recovery=group['PGM Recovery [%]'].min()
    max_recovery=group['PGM Recovery [%]'].max()
    min_flow=group['Q_org [L/time]'].min()
    max_flow=group['Q_org [L/time]'].max()
    print('Stages:', stages)
    print('Min Recovery:', min_recovery)
    print('Max Recovery:', max_recovery)
    print('Min Organic Flow Rate:', min_flow)
    print('Max Organic Flow Rate:', max_flow)
    stage_list.append(stages)
    min_recovery_list.append(min_recovery)
    max_recovery_list.append(max_recovery)
    min_flow_list.append(min_flow)
    max_flow_list.append(max_flow)
    # count+=1
    # print(count)
# print(df_random_filtered_sorted)



max_recov_list=df_constrained_nonrandom['Rh Recovery [%]'].tolist()
vol_flow_list=df_constrained_nonrandom['Q_org [L/time]'].tolist()
n_stages_list=df_constrained_nonrandom['Stages'].tolist()
desired_purity_Rh=df_constrained_nonrandom['Rh Relative Purity [%]'][0]

max_recov_arr=np.array(max_recov_list)
vol_flow_arr=np.array(vol_flow_list)
n_stages_arr=np.array(n_stages_list)
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
# for each of the stage numbers, I want to plot a vertical line that goes from the min recovery to the max recovery
for i in range(len(stage_list)):
    stages=stage_list[i]
    min_recovery=min_recovery_list[i]
    max_recovery=max_recovery_list[i]
    min_flow=min_flow_list[i]
    max_flow=max_flow_list[i]

    # plot vertical line for recovery range
    ax1.vlines(x=stages, ymin=min_recovery, ymax=max_recovery, color='tab:blue', alpha=0.5)
    # at the min and max, add a horizontal line like a cap to the vertical line to make it look like an error bar
    ax1.hlines(y=min_recovery, xmin=stages-0.2, xmax=stages+0.2, color='tab:blue', alpha=0.5)
    ax1.hlines(y=max_recovery, xmin=stages-0.2, xmax=stages+0.2, color='tab:blue', alpha=0.5)
    
    # plot vertical line for organic flow rate range
    ax2.vlines(x=stages, ymin=min_flow, ymax=max_flow, color='tab:red', alpha=0.5)
    # at the min and max, add a horizontal line like a cap to the vertical line to make it look like an error bar
    ax2.hlines(y=min_flow, xmin=stages-0.2, xmax=stages+0.2, color='tab:red', alpha=0.5)
    ax2.hlines(y=max_flow, xmin=stages-0.2, xmax=stages+0.2, color='tab:red', alpha=0.5)
# ax1.plot([], [], color='tab:blue', alpha=0.5, label='Recovery Uncertainty Range')
ax1.legend(lines, labels, loc='center right')
# add the vertical lines to the legend
ax1.plot([], [], color='tab:blue', alpha=0.5, label='Recovery Uncertainty Range')
ax1.legend(bbox_to_anchor = (1, 0.65))
ax2.plot([], [], color='tab:red', alpha=0.5, label='Organic Flow Uncertainty Range')

ax2.legend(bbox_to_anchor = (1, 0.5))

plt.title('Recovery at '+str(desired_purity_Rh)+'% Purity Rh and Organic Flow Vs Stage Number ddFc Ligand')

plt.tight_layout()
plt.show()

