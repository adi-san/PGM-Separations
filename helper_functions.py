import numpy as np
from scipy.optimize import fsolve
from numpy import ndarray
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, minimize, Bounds, LinearConstraint, NonlinearConstraint, curve_fit
import pandas as pd
import openpyxl
from scipy.stats import t
from datetime import datetime


# data extraction and fitting functions
def unary_langmuir_func(x, q_max, K):
    return q_max*K*x/(1+K*x)
def d_langmuir_func_d_q_max(x, q_max, K):
    return K*x/(1+K*x)
def d_langmuir_func_d_K(x, q_max, K):
    return q_max*x/(1+K*x)**2
def Jac_for_data_fitting(x_data, param_arr):
    q_max=param_arr[0]
    K=param_arr[1]
    Jac=np.ones((len(x_data),len(param_arr)))
    Jac[:,0]=d_langmuir_func_d_q_max(x_data, q_max, K)
    Jac[:,1]=d_langmuir_func_d_K(x_data, q_max, K)
    return Jac
def unary_langmuir_func_resid(theta_arr,x_data,y_data):
    q_max=theta_arr[0]
    K=theta_arr[1]
    return y_data-unary_langmuir_func(x_data, q_max, K)
# end data extraction and fitting functions


# from recovery_purity_tradeoff_analysis_countercurrent import C_in
# I want to solve for outlet conc
def denom_fcn(K_eq_arr: ndarray, C_arr: ndarray):
  return 1+np.dot(K_eq_arr,C_arr)
def q_func_langmuir_relation(q_max_arr: ndarray, K_eq_arr: ndarray, C_arr: ndarray):
  denom=denom_fcn(K_eq_arr,C_arr)
  q_arr=q_max_arr*K_eq_arr*C_arr/denom
  return q_arr
def ppm_mass_tot_to_M_concs(ppm_mass_tot: float ,rel_mol_frac_arr: ndarray,MW_arr: ndarray):
  # use basis of 1 mol (i.e. rel mol frac --> x_i mol for calculation purposes)
  dummy_mass_arr=rel_mol_frac_arr*MW_arr
  dummy_mass_tot=np.sum(dummy_mass_arr)
  rel_mass_frac_arr=dummy_mass_arr/dummy_mass_tot
  rel_mass_conc_arr= rel_mass_frac_arr*ppm_mass_tot/1000 #g/L
  M_concs=rel_mass_conc_arr/MW_arr
  return M_concs
# print(q_func_langmuir_relation(q_max_arr,K_Eq_arr,a))
# C_out is what we want to solve for. C_lig, Q_aq, Q_org are scalars that do not change
def mat_bal_func_stage(C_out: ndarray, C_lig: float, Q_aq: float, Q_org: float, q_in: ndarray,C_in: ndarray,q_max: ndarray,K_eq: ndarray):
  # wnat func_value to be 0
  func_value=Q_aq*C_in+Q_org*C_lig*q_in-Q_aq*C_out-Q_org*C_lig*q_func_langmuir_relation(q_max,K_eq,C_out)
  return func_value

def crosscurrent_model_fsolve(C_in, q_in, C_lig, Q_aq, Q_org, n_stages,q_max_arr,K_Eq_arr):
  q_arr=np.zeros((n_stages,len(C_in)))
  C_arr=np.zeros((n_stages+1,len(C_in)))
  C_arr[0,:]=C_in
  for i in range(n_stages):
    result=fsolve(mat_bal_func_stage, C_in, args=(C_lig,Q_aq,Q_org/n_stages,q_in,C_in,q_max_arr,K_Eq_arr),xtol=1e-6)
    C=result
    # C=fsolve(mat_bal_func_stage, C_in, args=(C_lig,Q_aq,Q_org/n_stages,q_in,C_in,q_max_arr,K_Eq_arr), xtol=1e-6)
    C_arr[i+1,:]=C
    # print()
    q_out=q_func_langmuir_relation(q_max_arr,K_Eq_arr,C) # Pass C.x instead of C
    C_in=C # C is already C.x from the line above
    q_arr[i,:]=q_out
  return C_arr, q_arr

def compute_recov_aq_arr(C_arr,n_stages):
  C_in=C_arr[0,:]
  recov_aq_arr=np.zeros((n_stages+1,len(C_in)))
# recov_org_arr=np.zeros((n_stages+1,len(C_in)))
  for i in range(n_stages+1):
    recov_aq_arr[i,:]=C_arr[i,:]/C_in*100
  return recov_aq_arr

def compute_Rh_purity_aq(C_arr):
  Rh_purity_aq=C_arr[-1,2]/np.sum(C_arr[-1,:])*100
  return Rh_purity_aq

def compute_Pd_purity_aq(C_arr):
  Pd_purity_aq=C_arr[-1,1]/np.sum(C_arr[-1,:])*100
  return Pd_purity_aq
def compute_Pt_purity_aq(C_arr):
  Pt_purity_aq=C_arr[-1,0]/np.sum(C_arr[-1,:])*100
  return Pt_purity_aq

def countercurrent_model(C_stages_1_to_n: ndarray, C_lig: float, Q_aq: float, Q_org: float,n_stages: int, q_in: ndarray,C_in: ndarray,q_max: ndarray,K_eq: ndarray):
  func_value=np.zeros(len(C_stages_1_to_n))

  # the below MUST be an integer
  n_comps=int(len(C_stages_1_to_n)/n_stages)
  scale_factor=10**0
  for i in range(n_stages):
    if i==0:
      func_value[0:n_comps]=(Q_aq*C_in+Q_org*C_lig*q_func_langmuir_relation(q_max, K_eq, C_stages_1_to_n[n_comps:2*n_comps])-Q_aq*C_stages_1_to_n[0:n_comps]-Q_org*C_lig*q_func_langmuir_relation(q_max, K_eq, C_stages_1_to_n[0:n_comps]))*scale_factor
      dummy_4_scale=func_value[n_comps-1].copy()
      func_value[n_comps-1]=dummy_4_scale/q_max[-1]
      # print(67)
    # else if we are on the last stage (recall that the indices are )
    elif i==n_stages-1:
      func_value[n_comps*(n_stages-1):n_comps*n_stages]=(Q_aq*C_stages_1_to_n[n_comps*(n_stages-2):n_comps*(n_stages-1)]
      +Q_org*C_lig*q_in
      -Q_aq*C_stages_1_to_n[n_comps*(n_stages-1):n_comps*n_stages]
      -Q_org*C_lig*q_func_langmuir_relation(q_max, K_eq, C_stages_1_to_n[n_comps*(n_stages-1):n_comps*(n_stages)]))*scale_factor
      dummy_4_scale=func_value[-1].copy()
      func_value[-1]=dummy_4_scale/q_max[-1]
    else:
      # suppose this is stage 2 for shits and giggles lol, so i=1
      # suppose this is stage 3 for shits and giggles lol, so i=2
      func_value[n_comps*(i):n_comps*(i+1)]=(Q_aq*C_stages_1_to_n[n_comps*(i-1):n_comps*(i)]
      +Q_org*C_lig*q_func_langmuir_relation(q_max, K_eq, C_stages_1_to_n[n_comps*(i+1):n_comps*(i+2)])
      -Q_aq*C_stages_1_to_n[n_comps*(i):n_comps*(i+1)]
      -Q_org*C_lig*q_func_langmuir_relation(q_max, K_eq, C_stages_1_to_n[n_comps*(i):n_comps*(i+1)]))*scale_factor
      dummy_4_scale=func_value[n_comps*(i+1)-1].copy()
      func_value[n_comps*(i+1)-1]=dummy_4_scale/q_max[-1]
  return func_value

# this is for the ddFc Ligand
def Rh_purity_resid_fcn_countercurrent(Q_org, 
                                       C_lig, 
                                       Q_aq,
                                       n_stages, 
                                       q_in,
                                       C_in,
                                       q_max,
                                       K_eq,
                                       purity_threshold):
  # print(C_in)
  low_b = np.zeros(len(C_in)*n_stages)
  up_b=np.ones(len(C_in)*n_stages)*np.inf
  opt_bounds = (low_b, up_b)
  C_stages_1_to_n_guess=np.tile(C_in,n_stages)
  C_arr=least_squares(countercurrent_model, C_stages_1_to_n_guess, args=(C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max,K_eq), ftol=1e-11, gtol=1e-11, xtol=1e-11,bounds=opt_bounds, method='trf')
  C_arr=C_arr.x
  if np.linalg.norm(countercurrent_model(C_arr, C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max,K_eq))>1e-6:
    print(f'Warning: The solution for Q_org={Q_org} may not have converged properly, as the l2 norm of the objective function is greater than 1e-6.')
  C_arr=C_arr.reshape(int(len(C_arr)/int(len(C_in))),int(len(C_in)))
  Rh_purity_aq=compute_Rh_purity_aq(C_arr)
  rel_resid=abs(Rh_purity_aq-purity_threshold)
  return rel_resid
# This is for the coeFc ligand
def Pt_purity_resid_fcn_countercurrent(Q_org, 
                                       C_lig, 
                                       Q_aq,
                                       n_stages, 
                                       q_in,
                                       C_in,
                                       q_max,
                                       K_eq,
                                       purity_threshold):
  # print(C_in)
  low_b = np.zeros(len(C_in)*n_stages)
  up_b=np.ones(len(C_in)*n_stages)*np.inf
  opt_bounds = (low_b, up_b)
  C_stages_1_to_n_guess=np.tile(C_in,n_stages)
  C_arr=least_squares(countercurrent_model, C_stages_1_to_n_guess, args=(C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max,K_eq), ftol=1e-11, gtol=1e-11, xtol=1e-11,bounds=opt_bounds, method='trf')
  C_arr=C_arr.x
  if np.linalg.norm(countercurrent_model(C_arr, C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max,K_eq))>1e-6:
    print(f'Warning: The solution for Q_org={Q_org} may not have converged properly, as the l2 norm of the objective function is greater than 1e-6.')
  C_arr=C_arr.reshape(int(len(C_arr)/int(len(C_in))),int(len(C_in)))
  Pt_purity_aq=compute_Pt_purity_aq(C_arr)
  rel_resid=abs(Pt_purity_aq-purity_threshold)
  return rel_resid

def water_flow_successive_stages(mol_uptake_arr_prev_stage: ndarray, 
                                 ppm_mass_tot_des: float, 
                                 MW_arr: ndarray, 
                                 org_flow_prev_stage: float,
                                 C_lig_prev_stage: float):
  rel_mol_fracs=mol_uptake_arr_prev_stage/np.sum(mol_uptake_arr_prev_stage)
  # dummy basis is 1 mol
  dummy_basis_mass=rel_mol_fracs*MW_arr
  rel_mass_fracs=dummy_basis_mass/np.sum(dummy_basis_mass)
  g_per_L_mass_concs=rel_mass_fracs*ppm_mass_tot_des/1000
  mol_concs=g_per_L_mass_concs/MW_arr
  mol_per_time_entering=mol_uptake_arr_prev_stage*C_lig_prev_stage*org_flow_prev_stage
  vol_flow_aq_in=mol_per_time_entering/mol_concs
  return vol_flow_aq_in[0]

# constrained purity analysis functions
def run_constrained_purity_analysis_countercurrent(C_in, 
                                                   q_in, 
                                                   C_lig, 
                                                   Q_aq,
                                                   Q_org_ig,  
                                                   q_max_arr, 
                                                   K_Eq_arr, 
                                                   desired_purity_PGM,
                                                   starting_stage,
                                                   highest_stage,
                                                   ligand,
                                                   PGM_labels,
                                                   MW_arr,
                                                   constrained_purity_func):
  n_stages_arr=np.linspace(starting_stage,highest_stage,highest_stage-starting_stage+1)
  n_stages_int_list=[]
  max_recov_list=[]
  vol_flow_list=[]
  C_out_list=[]
  q_out_list=[]
  for n_stages_float in n_stages_arr:
    n_stages = int(n_stages_float)
    n_stages_int_list.append(n_stages)
    low_b = np.zeros(len(C_in)*n_stages)
    up_b=np.ones(len(C_in)*n_stages)*np.inf
    opt_bounds = (low_b, up_b)

    result=fsolve(constrained_purity_func,Q_org_ig,args=(C_lig, Q_aq,n_stages, q_in,C_in,q_max_arr,K_Eq_arr,desired_purity_PGM), xtol=1e-6)
    
    Q_org_det=result[0]
  #   print(result)
    #test to see if the print is making runtime slower
    print(f'For n_stages={n_stages}, the required organic flow rate to achieve {desired_purity_PGM}% purity of Rh in the aqueous phase is {Q_org_det} L/time.')
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

    # append the last stage exit conc vec to c_out_list as we will make a DF of these and plot 
    # them to see how the exit concentrations change with number of stages at the point where we achieve the desired purity value
    C_out_list.append(C_counter_plot[-1,:])
    q_out=q_func_langmuir_relation(q_max_arr, K_Eq_arr, C_counter_plot[1,:])
    q_out_list.append(q_out)
    # price_Rh=price_arr[-1]*C_counter_plot[-1,-1] #$Rh/L feed
    # max_sales_Rh_list.append(price_Rh)
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

  Conc_out_labels=['C_'+label+'_out [M]' for label in PGM_labels]
  Conc_in_labels=['C_'+label+'_in [M]' for label in PGM_labels]
  q_in_labels=['q_'+label+'_in [mol PGM/mol ligand]' for label in PGM_labels]
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
  ppm_in=C_in*MW_arr*1000
  result['Q_org [L/time]']=vol_flow_arr
  result['PGM Recovery [%]']=max_recov_arr
  result['PPM Mass Total Feed [mg/L]']=np.sum(ppm_in)
  result['PGM Relative Purity [%]']=desired_purity_PGM
  result['Q_aq [L/time]']=Q_aq
  result['C_lig [L/time]']=C_lig
  result['Q_aq/(Q_org*C_lig) [1/M]']=Q_aq/(vol_flow_arr*C_lig)
  result['Ligand']=ligand
  # result['Max Sales Rh [$Rh/L feed basis]']=max_sales_Rh_list
  # result=result[['Stages']+Conc_in_labels+q_in_labels+Conc_out_labels+q_out_labels]
  # result=result[['Stages']+Conc_in_labels+q_in_labels+Conc_out_labels+q
  # print(result)
  return result