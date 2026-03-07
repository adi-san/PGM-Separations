from helper_functions import *



# Load the entire Excel file
all_sheets = pd.read_excel('unary_isothern_extraction_data_automation.xlsx', sheet_name=None)
confidence_level=0.95

with pd.ExcelFile('unary_isothern_extraction_data_automation.xlsx') as xlsx:
    # Get list of all sheet names
    names = xlsx.sheet_names
    # print(f"Sheets found: {names}")
    df_list=[]
    for sheet in names:
        df=pd.read_excel(xlsx, sheet)
        x_data=df['Ce (M)'].to_numpy()
        y_data=df['q (mol/mol)'].to_numpy()
        popt= least_squares(unary_langmuir_func_resid,x0=[0.01,2000],args=(x_data, y_data) ,ftol=1e-11, gtol=1e-11, xtol=1e-11)
        param_arr=popt.x
        y_hat=unary_langmuir_func(x_data,param_arr[0],param_arr[1])
        resid=y_hat-y_data
        # now I want to calculate the standard errors and confidence intervals for the fitted parameters using the Jacobian of the residuals at the solution, which is provided by the least_squares output
        Jac=Jac_for_data_fitting(x_data, param_arr)
        J_T_J=np.transpose(Jac) @ Jac
        J_T_J_inv=np.linalg.inv(J_T_J)
        n=len(x_data)
        p=len(param_arr)
        # we subtract another DOF to add more uncertainty
        DOF=n-(p+1)
        MSE=np.dot(resid,resid)/(DOF)
        covar_matrix=MSE*J_T_J_inv
        # collect the standard errors of the parameters, which are the square roots of the diagonal elements of the covariance matrix
        std_errors=np.sqrt(np.diag(covar_matrix))
        MOEs=std_errors*t.ppf((1 + confidence_level) / 2, DOF)
        # create a DF to store the results for this sheet, and then we will concatenate these DFs for all sheets to make a summary DF
        results_df=pd.DataFrame({'Parameter': ['q_max [mol PGM/mol Ligand]', 'K_eq [~]'], 'Value': param_arr, 'Standard Error': std_errors, 'Margin of Error '+str(confidence_level*100)+'% Confidence': MOEs})
        df_list.append(results_df)
        # print(results_df)
    # concatenate the results DFs for all sheets into a summary DF
    # I want the PGM_ligand to be the index of the summary DF, so I will set the index of each results_df to be the sheet name before concatenating
        results_df.set_index(pd.Index([sheet]*len(results_df)), inplace=True)
        # I want the index column to be named 'PGM_Ligand'
        results_df.index.name='PGM Ligand Pair'
    summary_df=pd.concat(df_list)
    print(summary_df)
    summary_df.to_csv('isotherm_parameters_summary.csv')



        # print(df)
        

    
    # # Read specific sheets
    # df1 = pd.read_excel(xls, 'Sheet1')
    # df2 = pd.read_excel(xls, 'Sheet2')
