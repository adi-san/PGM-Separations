from helper_functions import *

df_constrained_nonrandom=pd.read_csv('Constrained Purity Analysis Pt Pd Rh/constrained_purity_analysis_Pt_Pd_Rh_with_regressed_parameters.csv')
df_constrained_random=pd.read_csv('Pt_Pd_Rh_S1_S2_uncertainty_results.csv')

# Sort by the first column (assumed to be 'S1 Stages')
df_constrained_nonrandom = df_constrained_nonrandom.sort_values(by=df_constrained_nonrandom.columns[0])
df_constrained_random = df_constrained_random.sort_values(by=df_constrained_random.columns[0])
# Drop columns with prefix 'S2'
df_constrained_nonrandom_S1 = df_constrained_nonrandom.loc[:, ~df_constrained_nonrandom.columns.str.startswith('S2')]
df_constrained_random_S1 = df_constrained_random.loc[:, ~df_constrained_random.columns.str.startswith('S2')]
# print(df_constrained_nonrandom_S1)
# drop duplicate rows of df_constrained_nonrandom_S1
df_constrained_nonrandom_S1 = df_constrained_nonrandom_S1.drop_duplicates()
# drop the S1 prefix from all of the columns
df_constrained_nonrandom_S1.columns = df_constrained_nonrandom_S1.columns.str.replace('S1 ', '')
df_constrained_random_S1.columns = df_constrained_random_S1.columns.str.replace('S1 ', '')
# reconfigure indices
df_constrained_nonrandom_S1.reset_index(drop=True, inplace=True)
# print(df_constrained_nonrandom_S1)
df_constrained_random_S1.reset_index(drop=True, inplace=True)
# print(df_constrained_random_S1)


df_random_S1_filtered=df_constrained_random_S1[(df_constrained_random_S1['Q_org [L/time]']>0) & (df_constrained_random_S1['Q_org [L/time]']<1)]
# get the K_eq and q_max values from the rows we dropped, and see if they are in a reasonable range
dropped_rows=df_constrained_random_S1[~((df_constrained_random_S1['Q_org [L/time]']>0) & (df_constrained_random_S1['Q_org [L/time]']<1))]
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
    df_random_S1_filtered=df_random_S1_filtered[~((df_random_S1_filtered['K_eq_Pt']==K_eq_Pt) & (df_random_S1_filtered['q_max_Pt']==q_max_Pt) & (df_random_S1_filtered['K_eq_Pd']==K_eq_Pd) & (df_random_S1_filtered['q_max_Pd']==q_max_Pd) & (df_random_S1_filtered['K_eq_Rh']==K_eq_Rh) & (df_random_S1_filtered['q_max_Rh']==q_max_Rh))]

# print(df_random_filtered)
# sort the df_random_filtered by the number of stages
# df_random_filtered_sorted=df_random_filtered.sort_values(by='Stages')
# for each chunk where each chunk corresponds to a different number of stages, I want the min and max recovery nad the min and max organic flow rate

# count=0
stage_list=[]
min_recovery_list=[]
max_recovery_list=[]
low_quartile_recovery_list=[]
up_quartile_recovery_list=[]
min_flow_list=[]
max_flow_list=[]
low_quartile_org_flow_list=[]
up_quartile_org_flow_list=[]
for stages, group in df_random_S1_filtered.groupby('Stages'):
    min_recovery=group['PGM Recovery [%]'].min()
    max_recovery=group['PGM Recovery [%]'].max()
    recov_list=group['PGM Recovery [%]'].to_list()
    low_quartile_recovery=np.percentile(recov_list,25)
    up_quartile_recovery=np.percentile(recov_list,75)
    min_flow=group['Q_org [L/time]'].min()
    max_flow=group['Q_org [L/time]'].max()
    org_flow_list=group['Q_org [L/time]'].to_list()
    low_quartile_flow=np.percentile(org_flow_list,25)
    up_quartile_flow=np.percentile(org_flow_list,75)
    print('Stages:', stages)
    print('Min Recovery:', min_recovery)
    print('Max Recovery:', max_recovery)
    print('Low Quartile Recovery:', low_quartile_recovery)
    print('Upper Quartile Recovery:', up_quartile_recovery)
    print('Min Organic Flow Rate:', min_flow)
    print('Max Organic Flow Rate:', max_flow)
    print('Low Quartile Org Flow:', low_quartile_flow)
    print('Upper Quartile Org Flow:', up_quartile_flow)
    stage_list.append(stages)
    min_recovery_list.append(min_recovery)
    max_recovery_list.append(max_recovery)
    low_quartile_recovery_list.append(low_quartile_recovery)
    up_quartile_recovery_list.append(up_quartile_recovery)
    min_flow_list.append(min_flow)
    max_flow_list.append(max_flow)
    low_quartile_org_flow_list.append(low_quartile_flow)
    up_quartile_org_flow_list.append(up_quartile_flow)
    # count+=1
    # print(count)
# print(df_random_filtered_sorted)



max_recov_list=df_constrained_nonrandom_S1['PGM Recovery [%]'].tolist()
vol_flow_list=df_constrained_nonrandom_S1['Q_org [L/time]'].tolist()
n_stages_list=df_constrained_nonrandom_S1['Stages'].tolist()
desired_purity_Rh=df_constrained_nonrandom_S1['PGM Relative Purity [%]'][0]

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
    low_quartile_recovery=low_quartile_recovery_list[i]
    up_quartile_recovery=up_quartile_recovery_list[i]

    min_flow=min_flow_list[i]
    max_flow=max_flow_list[i]
    low_quartile_flow=low_quartile_org_flow_list[i]
    up_quartile_flow=up_quartile_org_flow_list[i]

    # plot vertical line for recovery range
    ax1.vlines(x=stages, ymin=min_recovery, ymax=max_recovery, color='tab:blue', alpha=0.5)
    # I want a vertical iqr at each point
    ax1.vlines(x=stages, ymin=low_quartile_recovery, ymax=up_quartile_recovery, color='tab:blue', alpha=0.2,linewidth=10)
    # at the min and max, add a horizontal line like a cap to the vertical line to make it look like an error bar
    ax1.hlines(y=min_recovery, xmin=stages-0.2, xmax=stages+0.2, color='tab:blue', alpha=0.5)
    ax1.hlines(y=max_recovery, xmin=stages-0.2, xmax=stages+0.2, color='tab:blue', alpha=0.5)
    
    # plot vertical line for organic flow rate range
    ax2.vlines(x=stages, ymin=min_flow, ymax=max_flow, color='tab:red', alpha=0.5)
    # I want a vertical iqr at each point
    ax2.vlines(x=stages, ymin=low_quartile_flow, ymax=up_quartile_flow, color='tab:red', alpha=0.2,linewidth=10)
    # at the min and max, add a horizontal line like a cap to the vertical line to make it look like an error bar
    ax2.hlines(y=min_flow, xmin=stages-0.2, xmax=stages+0.2, color='tab:red', alpha=0.5)
    ax2.hlines(y=max_flow, xmin=stages-0.2, xmax=stages+0.2, color='tab:red', alpha=0.5)
    # I 
# ax1.plot([], [], color='tab:blue', alpha=0.5, label='Recovery Uncertainty Range')
ax1.legend(lines, labels, loc='center right')
# add the vertical lines to the legend
ax1.plot([], [], color='tab:blue', alpha=0.5, label='Recovery Uncertainty Range')
ax1.legend(bbox_to_anchor = (1, 0.65))
ax2.plot([], [], color='tab:red', alpha=0.5, label='Organic Flow Uncertainty Range')

ax2.legend(bbox_to_anchor = (1, 0.5))

plt.title('Recovery at '+str(desired_purity_Rh)+'% Purity Rh and Organic Flow Vs Stage Number ddFc Ligand')

plt.tight_layout()
# save this figure as a png file
plt.savefig('Constrained Purity Analysis Pt Pd Rh/Pt_Pd_Rh_S1_uncertainty_analysis.png')
# plt.show()
# now I want to extend this analysis to the S2 data
df_random_filtered=df_constrained_random[(df_constrained_random['S1 Q_org [L/time]']>0) 
                                         & (df_constrained_random['S1 Q_org [L/time]']<1) 
                                         &(df_constrained_random['S2 Q_org [L/time]']>0) 
                                         & (df_constrained_random['S2 Q_org [L/time]']<1)]
# get the K_eq and q_max values from the rows we dropped, and see if they are in a reasonable range
dropped_rows=df_constrained_random[~((df_constrained_random['S1 Q_org [L/time]']>0) 
                                   & (df_constrained_random['S1 Q_org [L/time]']<1))
                                   & ~((df_constrained_random['S2 Q_org [L/time]']>0) 
                                       & (df_constrained_random['S2 Q_org [L/time]']<1))]
print('Dropped rows with invalid Q_org values:')
print(dropped_rows[['S2 K_eq_Pt', 'S2 q_max_Pt', 'S2 K_eq_Pd', 'S2 q_max_Pd', 'S2 K_eq_Rh', 'S2 q_max_Rh']])
# in df_random_filtered, i want to drop the rows containing any of the combos of K_eq and q_max values that are in the dropped_rows, since those combos seem to lead to invalid results
for index, row in dropped_rows.iterrows():
    K_eq_Pt=row['S2 K_eq_Pt']
    q_max_Pt=row['S2 q_max_Pt']
    K_eq_Pd=row['S2 K_eq_Pd']
    q_max_Pd=row['S2 q_max_Pd']
    K_eq_Rh=row['S2 K_eq_Rh']
    q_max_Rh=row['S2 q_max_Rh']
    df_random_filtered=df_random_filtered[~((df_random_filtered['S2 K_eq_Pt']==K_eq_Pt) & 
                                            (df_random_filtered['S2 q_max_Pt']==q_max_Pt) & 
                                            (df_random_filtered['S2 K_eq_Pd']==K_eq_Pd) & 
                                            (df_random_filtered['S2 q_max_Pd']==q_max_Pd) & 
                                            (df_random_filtered['S2 K_eq_Rh']==K_eq_Rh) & 
                                            (df_random_filtered['S2 q_max_Rh']==q_max_Rh))]
    K_eq_Pt=row['S1 K_eq_Pt']
    q_max_Pt=row['S1 q_max_Pt']
    K_eq_Pd=row['S1 K_eq_Pd']
    q_max_Pd=row['S1 q_max_Pd']
    K_eq_Rh=row['S1 K_eq_Rh']
    q_max_Rh=row['S1 q_max_Rh']
    df_random_filtered=df_random_filtered[~((df_random_filtered['S1 K_eq_Pt']==K_eq_Pt) & 
                                            (df_random_filtered['S1 q_max_Pt']==q_max_Pt) & 
                                            (df_random_filtered['S1 K_eq_Pd']==K_eq_Pd) & 
                                            (df_random_filtered['S1 q_max_Pd']==q_max_Pd) & 
                                            (df_random_filtered['S1 K_eq_Rh']==K_eq_Rh) & 
                                            (df_random_filtered['S1 q_max_Rh']==q_max_Rh))]
for S1_stages, group in df_random_filtered.groupby('S1 Stages'):
    # print('S1 Stages:', S1_stages)
    group.reset_index(drop=True, inplace=True)
    stage_list=[]
    min_recovery_list=[]
    max_recovery_list=[]
    low_quartile_recovery_list=[]
    up_quartile_recovery_list=[]
    min_flow_list=[]
    max_flow_list=[]
    low_quartile_org_flow_list=[]
    up_quartile_org_flow_list=[]
    for S2_stages, group2 in group.groupby('S2 Stages'):
        # print('S2 Stages:', S2_stages)
        group2.reset_index(drop=True, inplace=True)
        # for each S2 stage, I want to get the min and max recovery of Pt, the OG dataframe calculated this incorrectly
        group2['S2 PGM Recovery [%]'] = group2['S2 C_Pt_out [M]']*group2['S2 Q_aq [L/time]']/(group2['S1 C_Pt_in [M]']*group2['S1 Q_aq [L/time]'])*100
        min_recovery=group2['S2 PGM Recovery [%]'].min()
        max_recovery=group2['S2 PGM Recovery [%]'].max()
        recov_list=group2['S2 PGM Recovery [%]'].to_list()
        low_quartile_recovery=np.percentile(recov_list,25)
        up_quartile_recovery=np.percentile(recov_list,75)
        min_flow=group2['S2 Q_org [L/time]'].min()
        max_flow=group2['S2 Q_org [L/time]'].max()
        org_flow_list=group2['S2 Q_org [L/time]'].to_list()
        low_quartile_flow=np.percentile(org_flow_list,25)
        up_quartile_flow=np.percentile(org_flow_list,75)
        # print('Stages:', stages)
        # print('Min Recovery:', min_recovery)
        # print('Max Recovery:', max_recovery)
        # print('Low Quartile Recovery:', low_quartile_recovery)
        # print('Upper Quartile Recovery:', up_quartile_recovery)
        # print('Min Organic Flow Rate:', min_flow)
        # print('Max Organic Flow Rate:', max_flow)
        # print('Low Quartile Org Flow:', low_quartile_flow)
        # print('Upper Quartile Org Flow:', up_quartile_flow)
        stage_list.append(S2_stages)
        min_recovery_list.append(min_recovery)
        max_recovery_list.append(max_recovery)
        low_quartile_recovery_list.append(low_quartile_recovery)
        up_quartile_recovery_list.append(up_quartile_recovery)
        min_flow_list.append(min_flow)
        max_flow_list.append(max_flow)
        low_quartile_org_flow_list.append(low_quartile_flow)
        up_quartile_org_flow_list.append(up_quartile_flow)
        # print(group2)
    # print(group)
    df_constrained_nonrandom_S2= df_constrained_nonrandom[(df_constrained_nonrandom['S1 Stages']==S1_stages)]
    print(df_constrained_nonrandom_S2)
    df_constrained_nonrandom_S2['S2 PGM Recovery [%]']=df_constrained_nonrandom_S2['S2 C_Pt_out [M]']*df_constrained_nonrandom_S2['S2 Q_aq [L/time]']/(df_constrained_nonrandom_S2['S1 C_Pt_in [M]']*df_constrained_nonrandom_S2['S1 Q_aq [L/time]'])*100
    # df_constrained_nonrandom_S2['S2 PGM Recovery [%]']
    max_recov_list=df_constrained_nonrandom_S2['S2 PGM Recovery [%]'].tolist()
    vol_flow_list=df_constrained_nonrandom_S2['S2 Q_org [L/time]'].tolist()
    n_stages_list=df_constrained_nonrandom_S2['S2 Stages'].tolist()
    # desired_purity_Pt=df_constrained_nonrandom_S2['S2 PGM Relative Purity [%]'][0]
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
        low_quartile_recovery=low_quartile_recovery_list[i]
        up_quartile_recovery=up_quartile_recovery_list[i]

        min_flow=min_flow_list[i]
        max_flow=max_flow_list[i]
        low_quartile_flow=low_quartile_org_flow_list[i]
        up_quartile_flow=up_quartile_org_flow_list[i]

        # plot vertical line for recovery range
        ax1.vlines(x=stages, ymin=min_recovery, ymax=max_recovery, color='tab:blue', alpha=0.5)
        # I want a vertical iqr at each point
        ax1.vlines(x=stages, ymin=low_quartile_recovery, ymax=up_quartile_recovery, color='tab:blue', alpha=0.2,linewidth=10)
        # at the min and max, add a horizontal line like a cap to the vertical line to make it look like an error bar
        ax1.hlines(y=min_recovery, xmin=stages-0.2, xmax=stages+0.2, color='tab:blue', alpha=0.5)
        ax1.hlines(y=max_recovery, xmin=stages-0.2, xmax=stages+0.2, color='tab:blue', alpha=0.5)
        
        # plot vertical line for organic flow rate range
        ax2.vlines(x=stages, ymin=min_flow, ymax=max_flow, color='tab:red', alpha=0.5)
        # I want a vertical iqr at each point
        ax2.vlines(x=stages, ymin=low_quartile_flow, ymax=up_quartile_flow, color='tab:red', alpha=0.2,linewidth=10)
        # at the min and max, add a horizontal line like a cap to the vertical line to make it look like an error bar
        ax2.hlines(y=min_flow, xmin=stages-0.2, xmax=stages+0.2, color='tab:red', alpha=0.5)
        ax2.hlines(y=max_flow, xmin=stages-0.2, xmax=stages+0.2, color='tab:red', alpha=0.5)
        # I 
    # ax1.plot([], [], color='tab:blue', alpha=0.5, label='Recovery Uncertainty Range')
    ax1.legend(lines, labels, loc='center right')
    # add the vertical lines to the legend
    ax1.plot([], [], color='tab:blue', alpha=0.5, label='Recovery Uncertainty Range')
    ax1.legend(bbox_to_anchor = (1, 0.65))
    ax2.plot([], [], color='tab:red', alpha=0.5, label='Organic Flow Uncertainty Range')

    ax2.legend(bbox_to_anchor = (1, 0.5))

    plt.title('Recovery at '+str(95)+'% Purity Pt Stage 1 Stages: '+str(S1_stages))

    plt.tight_layout()
    # save this figure as a png file
    plt.savefig(f'Constrained Purity Analysis Pt Pd Rh/Uncertainty Plots/Pt_Pd_Rh_S1_{S1_stages}_stages_uncertainty_analysis_recovery.png')