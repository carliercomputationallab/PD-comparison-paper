# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:28:06 2023

@author: P70073624
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def calculate_aic(obj_fn, num_parameters):
    """
    Calculate the Akaike Information Criterion (AIC).

    Parameters:
    - obj_fn (numpy array or list): Array of objective function values.
    - num_parameters (int): Number of parameters in the model.

    Returns:
    - float: AIC value.
    """
    num_observation = 8
    likelihood = np.exp(-0.5 * obj_fn)  # Assuming obj_fn is the negative log-likelihood

    aic = 2 * num_parameters - 2 * np.log(likelihood)

    return aic.mean(), aic.std()

def calculate_aicc(obj_fn, num_parameters):
    """
    Calculate the corrected Akaike Information Criterion (AICc).

    Parameters:
    - obj_fn (numpy array or list): Array of objective function values.
    - num_parameters (int): Number of parameters in the model.
    - num_observations (int): Number of observations or sample size.

    Returns:
    - float: AICc value.
    """
    num_observations = 8
    likelihood = np.exp(-0.5 * obj_fn)  # Assuming obj_fn is the negative log-likelihood
    
    aic = 2 * num_parameters - 2 * np.log(likelihood)

    aicc = aic + (2 * num_parameters * (num_parameters + 1)) / (num_observations - num_parameters - 1)
    
    
    return aicc.mean(), aicc.std()

# Read the DataFrame from a file (adjust the file path and format accordingly)
df_TPM = pd.read_excel('fittedPS/fittedPS_TPM.xlsx', index_col = 0)
mean_L = df_TPM['L'].mean()
std_L = df_TPM['L'].std()
# Replace zero values in 'obj_fn' column with NaN
df_TPM['obj_fn'] = df_TPM['obj_fn'].replace(0, float('nan'))


df_Garred = pd.read_excel('fittedPS/fittedPS_GM.xlsx', index_col = 0)

# Replace all missing values in the entire DataFrame with NaN
df_Garred.replace(to_replace='missing_value', value=float('nan'), inplace=True)

df_W = pd.read_excel('fittedPS/fittedPS_WM_F0p5.xlsx', index_col = 0)
# df_W0 = pd.read_excel('fittedPS/fittedPS_WM_F0.xlsx', index_col = 0)
# df_W1 = pd.read_excel('fittedPS/fittedPS_WM_F1.xlsx', index_col = 0)

df_Graff = pd.read_excel('fittedPS/fittedPS_UGM.xlsx', index_col = 0)

df_Graff18 = pd.read_excel('fittedPS/fittedPS_UGM18.xlsx', index_col = 0)

df_simpleW = pd.read_excel('fittedPS/fittedPS_SWM.xlsx', index_col = 0)


df_list = [df_TPM, df_Garred, df_W, df_simpleW, df_Graff, df_Graff18]
model_names = ['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM-18']
###############################################
# PLOT OBJECTIVE FN AND COMP TIME
###############################################

# Initialize lists to store mean and std values
mean_obj_fn = []
std_obj_fn = []
mean_comp_time = []
std_comp_time = []

terms_to_remove = ['P1.csv', 'P14.csv', 'P15.csv', 'P18.csv', 'P19.csv', 'P22.csv', 'P23.csv', 'P28.csv', 'P29.csv'] #extraneal dwells we are focusing on physioneal
for i, df in enumerate(df_list):
    print(i)
    df = df[~df.iloc[:,0].isin(terms_to_remove)]
    df_list[i] = df

# Compute mean and std for each dataframe
for df in df_list:
    mean_obj_fn.append(df['obj_fn'].mean())
    std_obj_fn.append(df['obj_fn'].std())
    mean_comp_time.append(df['comp time'].mean())
    std_comp_time.append(df['comp time'].std())

# Create a figure
fig, ax = plt.subplots(1, 2, figsize=(8,4))
ax1 = ax[0]
ax2 = ax[1]
# Plot mean and std of 'obj_fn' on the left y-axis
x_values = np.arange(len(df_list))
ax1.errorbar(x_values, mean_obj_fn, yerr=std_obj_fn,  marker='o', color='k', capsize=5, ls = '')
# mean_W0 = df_W0['obj_fn'].mean()
# std_W0 = df_W0['obj_fn'].std()
# ax1.errorbar(1.8, mean_W0, yerr= std_W0, marker='d', color='r', capsize=5, ls = '' )
# mean_W1 = df_W1['obj_fn'].mean()
# std_W1 = df_W1['obj_fn'].std()
# ax1.errorbar(2.2, mean_W1, yerr= std_W1, marker='*', color='g', capsize=5, ls = '' )
ax1.set_ylabel('RMSE', color='k')
ax1.set_xticks(x_values)
ax1.set_xticklabels(model_names, rotation = 45)
ax1.text(0, 0.9, '(a)', transform = ax1.transAxes)
ax2.errorbar(x_values,np.array(mean_comp_time) / 1000, yerr=np.array(std_comp_time) / 1000, marker='o', color='k', capsize=5, ls = '')
ax2.set_ylabel('Computational time, *1000 s')
ax2.set_xticks(x_values)
ax2.set_xticklabels(model_names, rotation = 45)
ax2.text(0, 0.9, '(b)', transform = ax2.transAxes)

# Display the plot
plt.tight_layout()
plt.show()

plt.savefig('final-plots/obj-fn-comp-time.png', dpi = 600)
###############################################
# PLOT AIC
###############################################
# num_parameters = [6, 0, 12, 0, 11, 18]
# mean_AIC = []
# std_AIC = []
# mean_AICC = []
# std_AICC = []
# for i, df in enumerate(df_list):
#     mean_AIC.append(calculate_aic(df['obj_fn'], num_parameters[i])[0])
#     mean_AICC.append(calculate_aicc(df['obj_fn'], num_parameters[i])[0])
#     std_AIC.append(calculate_aic(df['obj_fn'], num_parameters[i])[1])
#     std_AICC.append(calculate_aicc(df['obj_fn'], num_parameters[i])[1])
    
# # Create a figure with two subplots
# fig, axs = plt.subplots(1, 2, figsize=(8,4), sharey=True)
# x_values = np.arange(len(df_list))
# # Subplot 1: AIC
# axs[0].errorbar(x_values, mean_AIC, yerr=std_AIC, label='AIC', marker='o', color='b', capsize=5, ls = '')
# axs[0].set_ylabel('AIC')
# axs[0].tick_params('y', colors='b')
# axs[0].set_title('Mean and Standard Deviation of AIC')
# axs[0].set_xticks(x_values)
# axs[0].set_xticklabels(model_names, rotation = 45)

# # Subplot 2: AICc
# axs[1].errorbar(x_values, mean_AICC, yerr=std_AICC, label='AICc', marker='o', color='r', capsize=5, ls = '')
# axs[1].set_xlabel('Dataframes')
# axs[1].set_ylabel('AICc')
# axs[1].tick_params('y', colors='r')
# axs[1].set_title('Mean and Standard Deviation of AICc')
# axs[1].set_xticks(x_values)
# axs[1].set_xticklabels(model_names, rotation = 45)

# # Display the plot
# plt.show()
# plt.savefig('final-plots/AIC.png', dpi = 600)    

###############################################
# PLOT MEAN PARAMETERS
###############################################
df_Garred = df_Garred.iloc[:,:-2].astype('float64')
df_Graff = df_Graff.iloc[:,:-2].astype('float64')
df_Graff18 = df_Graff18.iloc[:,:-2].astype('float64')
df_simpleW = df_simpleW.iloc[:,:-2].astype('float64')
df_TPM = df_TPM.iloc[:,:-6].astype('float64')
df_W = df_W.iloc[:,:-2].astype('float64')
df_list = [df_Graff18, df_TPM, df_Garred, df_W, df_simpleW, df_Graff]
terms_to_remove = ['P14.csv', 'P15.csv', 'P18.csv', 'P19.csv', 'P22.csv', 'P23.csv'] #extraneal dwells we are focusing on physioneal
for i, df in enumerate(df_list):
    print(i)
    df_list[i] = df[~df.index.isin(terms_to_remove)]
model_names = ['UGM-18','TPM', 'GM', 'WM', 'SWM', 'UGM']
result_df = pd.DataFrame(columns = model_names)
#%%
for i, df in enumerate(df_list):

    means = df.mean()
    
    stds = df.std()
    
    # Combine means and stds into a single dataframe
    result = pd.concat([means, stds], axis=1)
    result.columns = ['Mean', 'Std Dev']
    print(result)
    result_df[model_names[i]] = result.apply(lambda row: f"{row['Mean']:.2f} ± {row['Std Dev']:.2f}", axis=1)
    

result_df.loc['fct_urea', 'WM'] = 0
result_df.loc['fct_crea', 'WM'] = 0
result_df.loc['fct_glu', 'WM'] = 0
result_df.loc['fct_sodium', 'WM'] = 0
result_df.loc['fct_phosphate', 'WM'] = 0
result_df.loc['fct_potassium', 'WM'] = 0
result_df.loc['fct_urea', 'SWM'] = 0.5
result_df.loc['fct_crea', 'SWM'] = 0.5
result_df.loc['fct_glu', 'SWM'] = 0.5
result_df.loc['fct_sodium', 'SWM'] = 0.5
result_df.loc['fct_phosphate', 'SWM'] = 0.5
result_df.loc['fct_potassium', 'SWM'] = 0.5
result_df.loc['fct_urea', 'UGM'] = 1
result_df.loc['fct_crea', 'UGM'] = 1
result_df.loc['fct_glu', 'UGM'] = 1
result_df.loc['fct_sodium', 'UGM'] = 0.94
result_df.loc['fct_phosphate', 'UGM'] = 1
result_df.loc['fct_potassium', 'UGM'] = 0.96
result_df.loc['sico_urea', 'TPM'] = 0.963
result_df.loc['sico_crea', 'TPM'] = 0.958
result_df.loc['sico_glu', 'TPM'] = 0.948
result_df.loc['sico_sodium', 'TPM'] = 0.967
result_df.loc['sico_potassium', 'TPM'] = 0.961
result_df.loc['sico_phosphate', 'TPM'] = 0.961
result_df.loc['L', 'TPM'] = f"{mean_L:.2f} ± {std_L:.2f}"
result_df.loc['L', 'UGM'] = 0.3
result_df.loc['L', 'UGM-18'] = 0.3
result_df.to_excel('population_average_pigs.xlsx')
    
   

