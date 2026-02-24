import numpy as np
from scipy.optimize import fsolve
from numpy import ndarray
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, minimize, Bounds, LinearConstraint, NonlinearConstraint
import pandas as pd

from recovery_purity_tradeoff_analysis_countercurrent import C_in
# I want to solve for outlet conc
def denom_fcn(K_eq_arr: ndarray, C_arr: ndarray):
  return 1+np.dot(K_eq_arr,C_arr)
def q_func_langmuir_relation(q_max_arr: ndarray, K_eq_arr: ndarray, C_arr: ndarray):
  denom=denom_fcn(K_eq_arr,C_arr)
  q_arr=q_max_arr*K_eq_arr*C_arr/denom
  return q_arr
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

def Rh_purity_resid_fcn_countercurrent(Q_org, C_lig, Q_aq,n_stages, q_in,C_in,q_max,K_eq,purity_threshold):
  low_b = np.zeros(len(C_in)*n_stages)
  up_b=np.ones(len(C_in)*n_stages)*np.inf
  opt_bounds = (low_b, up_b)
  C_stages_1_to_n_guess=np.tile(C_in,n_stages)
  C_arr=least_squares(countercurrent_model, C_stages_1_to_n_guess, args=(C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max,K_eq), ftol=1e-11, gtol=1e-11, xtol=1e-11,bounds=opt_bounds, method='trf')
  C_arr=C_arr.x
  if np.linalg.norm(countercurrent_model(C_arr, C_lig, Q_aq, Q_org,n_stages, q_in,C_in,q_max,K_eq))>1e-8:
    print(f'Warning: The solution for Q_org={Q_org} may not have converged properly, as the l2 norm of the objective function is greater than 1e-8.')
  C_arr=C_arr.reshape(int(len(C_arr)/int(len(C_in))),int(len(C_in)))
  Rh_purity_aq=compute_Rh_purity_aq(C_arr)
  rel_resid=abs(Rh_purity_aq-purity_threshold)/purity_threshold
  return rel_resid