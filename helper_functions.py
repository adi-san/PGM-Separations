import numpy as np
from scipy.optimize import fsolve
from numpy import ndarray
import matplotlib.pyplot as plt
from scipy.optimize import minimize, Bounds

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