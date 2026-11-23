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
  ier_list=[]
  for n_stages_float in n_stages_arr:
    n_stages = int(n_stages_float)
    n_stages_int_list.append(n_stages)
    low_b = np.zeros(len(C_in)*n_stages)
    up_b=np.ones(len(C_in)*n_stages)*np.inf
    opt_bounds = (low_b, up_b)

    root, info, ier, mesg=fsolve(constrained_purity_func,Q_org_ig,args=(C_lig, Q_aq,n_stages, q_in,C_in,q_max_arr,K_Eq_arr,desired_purity_PGM), xtol=1e-6,full_output=True)
    
    Q_org_det=root[0]
  #   print(result)
    #test to see if the print is making runtime slower
    print(f'For n_stages={n_stages}, the required organic flow rate to achieve {desired_purity_PGM}% purity of Rh in the aqueous phase is {Q_org_det} L/time.')
    if ier!=1:
      print(f'Desired HK purity likely not met: ier= {ier}')
    # print(f'ier: {ier}')
    C_counter_ig=np.tile(C_in,n_stages)
    C_countercurrent=least_squares(countercurrent_model, C_counter_ig, args=(C_lig, Q_aq, Q_org_det,n_stages, q_in,C_in,q_max_arr,K_Eq_arr), ftol=1e-11, gtol=1e-11, xtol=1e-11,bounds=opt_bounds, method='trf')
    C_countercurrent_concs=C_countercurrent.x
  #   check the l2 norm of the objective function at the solution, which should be close to zero if the solution is good as mass is conserved across stages and thus the mass balance equations should be satisfied at the solution
    if np.linalg.norm(countercurrent_model(C_countercurrent_concs, C_lig, Q_aq, Q_org_det,n_stages, q_in,C_in,q_max_arr,K_Eq_arr))>1e-7:
      print(f'Warning: The solution for n_stages={n_stages} and Q_org={Q_org_det} may not have converged properly, as the l2 norm of the objective function is greater than 1e-7.')
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
    ier_list.append(ier)
    # price_Rh=price_arr[-1]*C_counter_plot[-1,-1] #$Rh/L feed
    # max_sales_Rh_list.append(price_Rh)
    # print(price_Rh)

    # below is an approach to initialize guesses: simple bisection method
    
    # if n_stages==2:
    #   Q_org_ig=Q_org_det/2
    # elif n_stages==3:
    #   Q_org_ig=Q_org_det/1.5
    # elif n_stages==4:
    #   Q_org_ig=Q_org_det*3/4
    # elif n_stages==5:
    #   Q_org_ig=Q_org_det*4/5
    Q_org_ig=Q_org_det*(n_stages-1)/n_stages
    
  
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
  result['ier']=ier_list
  # result['Max Sales Rh [$Rh/L feed basis]']=max_sales_Rh_list
  # result=result[['Stages']+Conc_in_labels+q_in_labels+Conc_out_labels+q_out_labels]
  # result=result[['Stages']+Conc_in_labels+q_in_labels+Conc_out_labels+q
  # print(result)
  return result

# ---------------------------------------------------------------------------
# Countercurrent model with incomplete ligand release and aqueous purge (recycle loop)
# ---------------------------------------------------------------------------
def recycle_q_in_factor(release_efficiency: float, purge_fraction: float):
  # q_in = (1-purge_fraction)*(1-release_efficiency)*q_out  (q_out = loading of organic leaving stage 1)
  return (1-purge_fraction)*(1-release_efficiency)

def countercurrent_with_recycle_resid(x: ndarray, C_lig: float, Q_aq: float, Q_org: float, n_stages: int,
                                      C_in: ndarray, q_max: ndarray, K_eq: ndarray, recycle_factor: float):
  # x = [C_stage1, ..., C_stage_n (n_stages*n_comps entries), q_in (n_comps entries)]
  n_comps=len(C_in)
  C_stages=x[:n_stages*n_comps]
  q_in=x[n_stages*n_comps:]
  resid_stages=countercurrent_model(C_stages, C_lig, Q_aq, Q_org, n_stages, q_in, C_in, q_max, K_eq)
  q_out=q_func_langmuir_relation(q_max, K_eq, C_stages[:n_comps])
  # recycle closure, scaled the same way as the Rh stage equations
  resid_recycle=(q_in-recycle_factor*q_out)/q_max[-1]
  return np.concatenate([resid_stages, resid_recycle])

def solve_countercurrent_with_recycle(Q_org, C_lig, Q_aq, n_stages, C_in, q_max, K_eq,
                                      release_efficiency, purge_fraction, guess=None, tol=1e-9):
  """Solves the stage concentrations AND the recycled ligand loading q_in simultaneously.
  Returns (C_mat [n_stages x n_comps], q_in [n_comps], l2 norm of residual)."""
  n_comps=len(C_in)
  rf=recycle_q_in_factor(release_efficiency, purge_fraction)
  n_unk=n_stages*n_comps+n_comps
  bounds=(np.zeros(n_unk), np.ones(n_unk)*np.inf)
  guesses=[]
  if guess is not None:
    guesses.append(guess)
  guesses.append(np.concatenate([np.tile(C_in,n_stages), np.zeros(n_comps)]))
  guesses.append(np.concatenate([np.tile(C_in*0.5,n_stages), 0.5*q_max*rf]))
  guesses.append(np.concatenate([np.tile(C_in*0.1,n_stages), 0.1*q_max*rf]))
  best=None
  for g in guesses:
    sol=least_squares(countercurrent_with_recycle_resid, g, args=(C_lig,Q_aq,Q_org,n_stages,C_in,q_max,K_eq,rf),
                      ftol=1e-12, gtol=1e-12, xtol=1e-12, bounds=bounds, method='trf')
    nrm=np.linalg.norm(sol.fun)
    if best is None or nrm<best[1]:
      best=(sol.x,nrm)
    if nrm<tol:
      break
  x,nrm=best
  return x[:n_stages*n_comps].reshape(n_stages,n_comps), x[n_stages*n_comps:], nrm, x

def constrained_purity_analysis_ddFc_with_release_acct(C_in, C_lig, Q_aq, n_stages, PGM_labels, q_max_arr, K_Eq_arr,
                                                       desired_purity_Rh, release_efficiency, purge_fraction,
                                                       price_arr, Q_org_ig=0.3, Q_org_min=0.02, Q_org_max=20, n_grid=40):
  """Finds the smallest organic flowrate that gives the desired Rh purity (%) in the aqueous exit when the
  organic phase is recycled with incomplete release and a purge. The recycled loading q_in is solved for, not specified.

  Purity is not monotone in Q_org and may never reach the target (the recycled q_in puts a floor on how clean the
  aqueous phase can get), so Q_org is found by scanning a log grid and refining the first crossing with brentq.
  Q_org_ig only orders the grid (it is kept for call compatibility with the previous scripts).

  Returns C_out, q_out, max_recov, vol_flow, max_sales_Rh, details (dict).
  If the target is unreachable, feasible=False and the Q_org with the highest purity on the grid is returned."""
  from scipy.optimize import brentq
  n_comps=len(C_in)
  state={'guess':None}
  def purity_at(Q_org):
    C_mat,q_in,nrm,x=solve_countercurrent_with_recycle(Q_org,C_lig,Q_aq,n_stages,C_in,q_max_arr,K_Eq_arr,
                                                       release_efficiency,purge_fraction,guess=state['guess'])
    if nrm<1e-7:
      state['guess']=x
    return compute_Rh_purity_aq(C_mat), nrm

  Q_grid=np.logspace(np.log10(Q_org_min),np.log10(Q_org_max),n_grid)
  purities=[]
  for Q in Q_grid:
    p,nrm=purity_at(Q)
    purities.append(p if nrm<1e-7 else np.nan)
  purities=np.array(purities)
  crossing=np.where(purities>=desired_purity_Rh)[0]
  feasible=len(crossing)>0
  if feasible:
    k=crossing[0]
    if k==0:
      Q_org_det=Q_grid[0]
    else:
      Q_org_det=brentq(lambda Q: purity_at(Q)[0]-desired_purity_Rh, Q_grid[k-1], Q_grid[k], xtol=1e-10, rtol=1e-10)
  else:
    Q_org_det=Q_grid[np.nanargmax(purities)]
  C_mat,q_in,nrm,x=solve_countercurrent_with_recycle(Q_org_det,C_lig,Q_aq,n_stages,C_in,q_max_arr,K_Eq_arr,
                                                     release_efficiency,purge_fraction,guess=state['guess'])
  if nrm>1e-7:
    print(f'Warning: n_stages={n_stages}, Q_org={Q_org_det}: residual norm {nrm:.2e} (may not have converged).')
  C_out=C_mat[-1,:]
  q_out=q_func_langmuir_relation(q_max_arr, K_Eq_arr, C_mat[0,:])
  max_recov=C_out[-1]/C_in[-1]*100
  max_sales_Rh=price_arr[-1]*C_out[-1] # $Rh/L feed
  purity=compute_Rh_purity_aq(C_mat)
  # overall balance with the recycle: feed = aqueous exit + (released + purged) PGM leaving the loop
  rf=recycle_q_in_factor(release_efficiency,purge_fraction)
  removed=Q_org_det*C_lig*(1-rf)*q_out
  mat_bal=Q_aq*C_in-Q_aq*C_out-removed
  details={'feasible':feasible,'purity':purity,'max_purity_on_grid':np.nanmax(purities),
           'C_stages':C_mat,'q_in':q_in,'q_out':q_out,'residual_norm':nrm,'removed_by_release_and_purge':removed,
           'overall_mat_bal':mat_bal,'Q_grid':Q_grid,'purity_grid':purities}
  return C_out, q_out, max_recov, Q_org_det, max_sales_Rh, details
