from helper_functions import *
# from run_countercurrent_simulations import K_Eq_arr, q_max_arr, C_in, q_in, C_lig, Q_aq, test_flowrates
# S2 of 

n_stages_arr=np.array([5])
q_in=np.array([0.0,0.0,0.0])
# total_conc_ppm_mass=500 # mg/L
# rel_mol_frac_arr=np.array([1/3, 1/3, 1/3]) #Pt, Pd, Rh
# MW_arr=np.array([195.084,192.22,102.91])
desired_purity_Pt=95

C_in=np.array([0.00127250930814378,0.00129714681171193,0.0000234790980741356])

Q_aq=0.773332585876166 # L/time
C_lig=0.1 #mols of ligand/L solution
test_flowrates=np.linspace(0.0,.25,50)

q_max_arr = np.array([0.272571356225773,0.409453841140154,0.0031092147446044]) #mol PGM/mol deFc
K_Eq_arr = np.array([349.540544019572,950.330400338052,369.020998529195]) # truly dimensionless parameters,

for n_stages in n_stages_arr:
    Pt_recov_list=[]
    Pt_purity_list=[]
    for Q_org in test_flowrates:
        C_arr_cross, q_arr_cross = crosscurrent_model_fsolve(C_in, q_in, C_lig, Q_aq, Q_org, n_stages,q_max_arr,K_Eq_arr)
        low_b = np.zeros(len(C_in)*n_stages)
        up_b=np.ones(len(C_in)*n_stages)*np.inf

        C_arr_cross_minus_in=C_arr_cross[1:,:]
        C_counter_ig=C_arr_cross_minus_in.flatten()

        opt_bounds = (low_b, up_b)
        C_countercurrent=least_squares(countercurrent_model, C_counter_ig, args=(C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max_arr,K_Eq_arr), ftol=1e-11, gtol=1e-11, xtol=1e-11,bounds=opt_bounds, method='trf')
        C_countercurrent_concs=C_countercurrent.x
        # check tle l@ norm of the objective function at the solution, which should be close to zero if the solution is good
        # print(f'For n_stages={n_stages} and Q_org={Q_org}, the l2 norm of the objective function at the solution is: {np.linalg.norm(countercurrent_model(C_countercurrent_concs, C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max_arr,K_Eq_arr))}')
        if np.linalg.norm(countercurrent_model(C_countercurrent_concs, C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max_arr,K_Eq_arr))>1e-7:
            print(f'Warning: The solution for n_stages={n_stages} and Q_org={Q_org} may not have converged properly, as the l2 norm of the objective function is greater than 1e-7.')
        C_countercurrent_mat=C_countercurrent_concs.reshape(int(len(C_countercurrent_concs)/int(len(C_in))),int(len(C_in)))

        C_counter_plot=np.vstack((C_in,C_countercurrent_mat)) # Fix: Use np.vstack for vertical concatenation

        recov_aq_arr=compute_recov_aq_arr(C_counter_plot,n_stages)
        Pt_purity_aq=compute_Pt_purity_aq(C_counter_plot)
        Pt_recov_list.append(recov_aq_arr[-1,0])
        Pt_purity_list.append(Pt_purity_aq)
    Rh_recov_arr=np.array(Pt_recov_list)
    Rh_purity_arr=np.array(Pt_purity_list)
    fig, ax1 = plt.subplots()

    # Recovery (left y-axis)
    ax1.plot(test_flowrates, Rh_recov_arr, color='tab:blue', label='Recovery')
    ax1.set_xlabel('Organic Flowrate (L/time)')
    ax1.set_ylabel('Recovery (%)', color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    # Purity (right y-axis)
    ax2 = ax1.twinx()
    ax2.plot(test_flowrates, Rh_purity_arr, color='tab:red', label='Purity')
    ax2.plot(test_flowrates, np.ones(len(test_flowrates))*desired_purity_Pt, color='tab:red',linestyle='--', label='Purity Threshold')
    ax2.set_ylabel('Pt Purity (%)', color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    # Optional: combined legend
    lines = ax1.get_lines() + ax2.get_lines()
    labels: list[str] = [str(line.get_label()) for line in lines]    
    ax1.legend(lines, labels, loc='center right')

    plt.title('Recovery and Purity Vs Org Flow at '+str(n_stages)+' Stages')

    plt.tight_layout()
    plt.show()