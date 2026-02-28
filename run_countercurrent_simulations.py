from helper_functions import *


col_labels = ['Pt', 'Pd', 'Rh']

q_in=np.array([0,0,0])
total_conc_ppm_mass=500 # mg/L
rel_mol_frac_arr=np.array([.45, .45, .1]) #Pt, Pd, Rh
MW_arr=np.array([195.084,106.42,102.91])
C_in=ppm_mass_tot_to_M_concs(total_conc_ppm_mass ,rel_mol_frac_arr,MW_arr)
n_stages=3
Q_aq=1 # L/time
C_lig=0.1 #mols of ligand/L solution
Q_org=0.7

eps=1e-11

q_max_arr = np.array([0.283186498156831,0.569456410590181,0.00376984198947194]) #mol PGM/mol ddFc
K_Eq_arr = np.array([1822.14319447879,2662.60543887326, 2401.64491456341]) # truly dimensionless parameters, I may need to scale these by powers of 10 if convergence proves to be tricky

# low_b=np.zeros(len(C_in)*n_stages)
low_b = np.zeros(len(C_in)*n_stages)
up_b=np.ones(len(C_in)*n_stages)*np.inf

C_arr_cross, q_arr_cross = crosscurrent_model_fsolve(C_in, q_in, C_lig, Q_aq, Q_org, n_stages,q_max_arr,K_Eq_arr)
C_arr_cross_minus_in=C_arr_cross[1:,:]
C_counter_ig=C_arr_cross_minus_in.flatten()

opt_bounds = (low_b, up_b)

C_countercurrent=least_squares(countercurrent_model, C_counter_ig, args=(C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max_arr,K_Eq_arr), ftol=eps, gtol=eps, xtol=eps,bounds=opt_bounds, method='trf')
print(C_countercurrent)
C_countercurrent_concs=C_countercurrent.x

C_arr_cross_df=pd.DataFrame(C_arr_cross,columns=col_labels)
C_arr_cross_df.plot(linestyle='--',marker='v')
plt.title('Crosscurrent Conc Vs Stages')
plt.xlabel("Stage Number")
plt.ylabel("Conc Profile [M]")
plt.show()

C_countercurrent_mat=C_countercurrent_concs.reshape(int(len(C_countercurrent_concs)/int(len(C_in))),int(len(C_in)))

C_counter_plot=np.vstack((C_in,C_countercurrent_mat)) # Fix: Use np.vstack for vertical concatenation
# print(C_countercurrent_mat)

print(f'Initial guess for countercurrent model: {C_counter_ig}')
print(f'Converged concentrations for countercurrent model: {C_countercurrent_concs}')

C_counter_plot_df=pd.DataFrame(C_counter_plot,columns=col_labels)

C_counter_plot_df.plot(linestyle='--',marker='v')
plt.title('Countercurrent Conc Vs Stages')
plt.xlabel("Stage Number")
plt.ylabel("Conc Profile [M]")
plt.show()

resid=countercurrent_model(C_countercurrent_concs, C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max_arr,K_Eq_arr)
print("L2 Norm of Residual", (np.linalg.norm(resid, ord=np.inf))**1)
print("Residuals", resid)
n_comps = len(C_in)

uptake_list=[q_func_langmuir_relation(q_max_arr, K_Eq_arr, C_countercurrent_concs[i*n_comps:(i+1)*n_comps]) for i in range(n_stages)]
uptake_mat=np.array(uptake_list).reshape(int(len(C_countercurrent_concs)/int(len(C_in))),int(len(C_in)))
uptake_mat_df=pd.DataFrame(uptake_mat, columns=col_labels)
uptake_mat_df.plot(linestyle='--',marker='v')
plt.title('Countercurrent Uptake Vs Stages')
plt.xlabel("Stage Number")
plt.ylabel("Uptake [mol PGM/mol ddFc]")
plt.show()

# print()
q_max_repeat=np.tile(q_max_arr,n_stages)
q_max_mat=np.array(q_max_repeat).reshape(int(len(C_countercurrent_concs)/int(len(C_in))),int(len(C_in)))
print(q_max_mat)
print("Purity of Rh in Aqueous Exit", C_countercurrent_concs[-1]/np.sum(C_countercurrent_concs[-3:]))
print("Recovery of Rh in Aqueous Exit", C_countercurrent_concs[-1]/C_in[-1])
mat_bal=C_countercurrent_mat[-1,:]*Q_aq+Q_org*C_lig*uptake_mat[0,:]
print('Material In [mol/hr] '+str(C_in*Q_aq))
print('Material Out [mol/hr] '+str(mat_bal))
frac_coverage_mat=uptake_mat/q_max_mat
new_col=[1-np.sum(frac_coverage_mat[i,:]) for i in range(n_stages)]
total_frac_coverage_mat=np.column_stack((frac_coverage_mat, new_col))

total_frac_coverage_mat_df=pd.DataFrame(total_frac_coverage_mat, columns=['Pt','Pd','Rh','Empty Sites'])
total_frac_coverage_mat_df.plot(linestyle='--',marker='v')
plt.title('Fractional Coverage Vs Stages')
plt.xlabel("Stage Number")
plt.ylabel("Fractional Coverage")
plt.show()
# print()