from helper_functions import *

# I want to generate random parameters for the isotherm parameters and see how the required organic flow rate changes to achieve a particular purity of Rh in the aqueous phase, as a function of the number of stages. I will use the same model and solver as in the constrained purity analysis, but I will loop over different sets of isotherm parameters.
df=pd.read_csv('isotherm_parameters_summary.csv')

# for each row, I want to generate 100 random samples of the parameters based on the standard errors
# the random samples will be drawn from a uniform distribution within the standard error for each parameter, which is a common way to represent uncertainty in parameters
num_samples=10000
random_params_list=[]
for index, row in df.iterrows():
    param_name=row['Parameter']
    param_value=row['Value']
    std_error=row['Standard Error']
    # MOE=row['Margin of Error 95.0% Confidence']
    # generate random samples for this parameter
    random_samples=np.random.uniform(param_value-std_error, param_value+std_error, num_samples)
    random_params_list.append(random_samples)
# now I have a list of arrays, where each array contains 100 random samples for a particular parameter. I want to combine these into a single array where each row is a set of parameters (q_max and K_eq for each PGM), and there are 100 rows corresponding to the 100 random samples.
random_params_array=np.array(random_params_list)
print(random_params_array)
# write a df with the random parameters, structured like the original df but with 100 rows for each parameter, and an additional column for the sample number
# i want to append the random parameters to the original df, so that I have a single df with all the original parameters and the random samples, which I can then use to run the countercurrent simulations for each set of parameters and analyze how the required organic flow rate changes with the isotherm parameters
# I will create new columns in the original df for each sample, and fill them with the random parameters
new_column_names=[]
for i in range(num_samples):
    new_column_names.append('Random Sample '+str(i+1))
# print(new_column_names)
df_rand=pd.DataFrame(random_params_array, columns=new_column_names)
df_combined=pd.concat([df, df_rand], axis=1)
print(df_combined)
df_combined.to_csv('isotherm_parameters_with_random_samples.csv', index=False)