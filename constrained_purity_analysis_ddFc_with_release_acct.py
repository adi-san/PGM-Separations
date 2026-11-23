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
# q_in=np.array([0,0,0])

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

"""
This script should be similar to the run_constrained_purity_analysis_countercurrent_Rh_From_Pt_Pd.py script, 
but it will account for the fact that the ligand is not fully released and that some of the aqueous phase is purged. 
The release efficiency and purge fraction are specified above. The rest of the code is similar to the previous script.
Instead of varying just the organic flowrate that gets us a 95% purity of Rh, we also need to vary the q_in values for Pt, Pd, and Rh
# we know that q_in_species_i= (1-purge_fraction)*(1-release_efficiency)*q_out_species_i, 
# so we can use this to calculate the q_in values for each species based on the q_out values from the previous stage.
# might need to write a function that allows for us to solve for the q_in values based on the q_out values as well as the purity metric.
"""

import os
plt.switch_backend('Agg') # save images without opening windows

out_dir='release_acct_stage_profiles'
os.makedirs(out_dir,exist_ok=True)

def plot_stage_profiles(details, n_stages, Q_org, release_eff, purge_frac, purity_target, ligand, PGM_labels, C_in, q_max_arr, out_dir):
    C_mat=details['C_stages']
    q_in=details['q_in']
    stages=np.arange(1,n_stages+1)
    q_mat=np.array([q_func_langmuir_relation(q_max_arr,K_Eq_arr,C_mat[i,:]) for i in range(n_stages)])
    frac_cov=q_mat/q_max_arr
    empty=1-np.sum(frac_cov,axis=1)
    fig,axs=plt.subplots(1,3,figsize=(17,4.8))
    for j,lab in enumerate(PGM_labels):
        axs[0].plot(np.arange(0,n_stages+1),np.vstack((C_in,C_mat))[:,j],'--v',label=lab)
        axs[1].plot(stages,q_mat[:,j],'--v',label=lab)
        axs[2].plot(stages,frac_cov[:,j],'--v',label=lab)
    axs[2].plot(stages,empty,'--v',color='k',label='Empty Sites')
    axs[0].set_title('Aqueous Concentration (stage 0 = feed)'); axs[0].set_ylabel('Conc [M]')
    axs[1].set_title('Ligand Uptake Leaving Each Stage'); axs[1].set_ylabel('Uptake [mol PGM/mol '+ligand+']')
    axs[2].set_title('Fractional Coverage'); axs[2].set_ylabel('Fractional Coverage')
    for ax in axs:
        ax.set_xlabel('Stage Number (aqueous enters stage 1, organic enters stage '+str(n_stages)+')')
        ax.legend(); ax.grid(alpha=0.3)
    qin_txt=', '.join(f'{lab}={v:.3g}' for lab,v in zip(PGM_labels,q_in))
    fig.suptitle(f'{ligand}, {n_stages} stages, Q_org={Q_org:.4g} L/time, Rh purity={details["purity"]:.2f}% (target {purity_target}%), '
                 f'release={release_eff}, purge={purge_frac}\nsolved recycled q_in [mol PGM/mol {ligand}]: {qin_txt}',fontsize=10)
    fig.tight_layout()
    fname=os.path.join(out_dir,f'stage_profiles_{ligand}_{n_stages}stages_release{release_eff}_purge{purge_frac}.png')
    fig.savefig(fname,dpi=150)
    plt.close(fig)
    return fname

# the first scenario is the one specified above; the second shows a release efficiency at which the target can actually be reached
scenarios=[release_efficiency, 0.9]
all_results=[]
for rel_eff in scenarios:
    for n_stages_float in n_stages_arr:
        n_stages_int=int(n_stages_float)
        C_out,q_out,max_recov,vol_flow,max_sales_Rh,details=constrained_purity_analysis_ddFc_with_release_acct(C_in,C_lig,Q_aq,n_stages_int,PGM_labels,q_max_arr,K_Eq_arr,desired_purity_Rh,rel_eff,purge_fraction,price_arr,Q_org_ig)
        if details['feasible']:
            print(f'release={rel_eff}, n_stages={n_stages_int}: Q_org={vol_flow:.5g} L/time gives {details["purity"]:.2f}% Rh purity, Rh recovery={max_recov:.2f}%, q_in={details["q_in"]}', flush=True)
            plot_stage_profiles(details,n_stages_int,vol_flow,rel_eff,purge_fraction,desired_purity_Rh,ligand,PGM_labels,C_in,q_max_arr,out_dir)
        else:
            print(f'release={rel_eff}, n_stages={n_stages_int}: {desired_purity_Rh}% Rh purity is NOT reachable (max purity on Q_org grid = {details["max_purity_on_grid"]:.2f}%)', flush=True)
        row={'Release Efficiency':rel_eff,'Purge Fraction':purge_fraction,'Stages':n_stages_int,'Feasible':details['feasible'],
             'Q_org [L/time]':vol_flow if details['feasible'] else np.nan,'Rh Purity [%]':details['purity'] if details['feasible'] else np.nan,
             'Max Rh Purity on Q_org grid [%]':details['max_purity_on_grid'],'Rh Recovery [%]':max_recov if details['feasible'] else np.nan,
             'Max Sales Rh [$Rh/L feed basis]':max_sales_Rh if details['feasible'] else np.nan}
        for j,lab in enumerate(PGM_labels):
            row['C_'+lab+'_out [M]']=C_out[j]
            row['q_'+lab+'_out [mol PGM/mol '+ligand+']']=q_out[j]
            row['q_'+lab+'_in (solved) [mol PGM/mol '+ligand+']']=details['q_in'][j]
        all_results.append(row)

results_df=pd.DataFrame(all_results)
results_df.to_csv(os.path.join(out_dir,'release_acct_results_'+ligand+'.csv'),index=False)

# purity vs Q_org for the specified scenario: shows why the target is (un)reachable
fig,ax=plt.subplots(figsize=(7,5))
for n_stages_int in [3,6,12]:
    if n_stages_int>highest_stage: continue
    _,_,_,_,_,d=constrained_purity_analysis_ddFc_with_release_acct(C_in,C_lig,Q_aq,n_stages_int,PGM_labels,q_max_arr,K_Eq_arr,desired_purity_Rh,release_efficiency,purge_fraction,price_arr,Q_org_ig)
    ax.semilogx(d['Q_grid'],d['purity_grid'],'-o',ms=3,label=f'{n_stages_int} stages')
ax.axhline(desired_purity_Rh,color='r',ls='--',label='target')
ax.set_xlabel('Q_org [L/time]'); ax.set_ylabel('Rh purity in aqueous exit [%]')
ax.set_title(f'Rh purity vs organic flow, release={release_efficiency}, purge={purge_fraction}'); ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(out_dir,f'purity_vs_Qorg_release{release_efficiency}_purge{purge_fraction}.png'),dpi=150); plt.close(fig)
