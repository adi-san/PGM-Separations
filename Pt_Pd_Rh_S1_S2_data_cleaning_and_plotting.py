from helper_functions import *
# Load the data
num_samples=500

# for i in range(num_samples):
#     # Load the data
#     df = pd.read_csv(f"Constrained Purity Analysis Pt Pd Rh/Uncertainty Results/constrained_purity_analysis_Pt_Pd_Rh_with_regressed_parameters_random_sample_{i+1}.csv")
    
#     # # Clean the data
#     print(df)
    
#     # # Plot the data
#     # plot_data(data, i)
df_list=[pd.read_csv(
    f"Constrained Purity Analysis Pt Pd Rh/Uncertainty Results/constrained_purity_analysis_Pt_Pd_Rh_with_regressed_parameters_random_sample_{i+1}.csv") 
    for i in range(num_samples)]
df=pd.concat(df_list, ignore_index=True)
print(df)
df.to_csv("Pt_Pd_Rh_S1_S2_uncertainty_results.csv", index=False)