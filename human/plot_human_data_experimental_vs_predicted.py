# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 09:52:51 2024

@author: P70073624
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

datafile = pd.read_excel('PD_in_silico_csv_export_/PD_in_silico_export_20240723.xlsx')
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)
patientlist = ['PDPC_059', 'PDPC_236','PDPC_268','PDPC_541','PDPC_572','PDPC_800','PDPC_835','PDPC_890','PDPC_946']


#Experimental
data = {}
data_V = {}

for i,ind in enumerate(datafile.index):
    df_cr = datafile.loc[ind,['Creatinine_t20', 'Creatinine_t120', 'Creatinine_t240']]*0.001
    df_cr.index = [20, 120, 240]
    
    df_glu = datafile.loc[ind,['Gluc_t20', 'Gluc_t120', 'Gluc_t240']]
    df_glu.index = [20, 120, 240]
    
    df_ur = datafile.loc[ind,['Ur_t20', 'Ur_t120', 'Ur_t240']]
    df_ur.index = [20, 120, 240]
    
    df_na = datafile.loc[ind,['Na20', 'Na120', 'Na240']]
    df_na.index = [20, 120, 240]
    
    df_cd = pd.concat([df_ur, df_cr, df_glu, df_na], axis = 1)
    df_cd.columns = ["Urea", "Creatinine", "Glucose", "Sodium"] #mmol/L
    
    #fill volume from datafile    
    Vfill = datafile.loc[ind, 'Vin_t0']
    
    # calculations of residual volume 

        
    RV0 = datafile.loc[ind, "RV0"]
    RV240 = datafile.loc[ind, "RV240"]
        
    print(RV0, RV240, flush= True)
    ur_0 = [datafile.loc[ind, 'Ur_night']*RV0/(Vfill+RV0) if datafile.loc[ind, 'Ur_night'] > 0 else 0]
    cr_0 = [datafile.loc[ind, 'Cr_night']*RV0/(Vfill+RV0)*0.001 if datafile.loc[ind, 'Cr_night'] > 0 else 0]
    glu_0 = [214 if datafile.loc[ind, 'Gluc_PET'] == 3.86 else  126 if datafile.loc[ind, 'Gluc_PET'] == 2.27 else  75 ]
    na_0 = [132]
    df_cd.loc[0] = np.concatenate((ur_0, cr_0, glu_0, na_0))
    
    # Volume from datafiles
    df_V = np.array(datafile.loc[ind,['Vin_t0', 'Vdrain_t240']]) #mL
    df_V[0] += RV0 #assuming a residual volume of 200 ml
    df_V[1] += RV240
    #Linear interpolation to find values of V at all times
    f_V = interp1d([0,240], df_V)
    interpolated_V = f_V(range(0,241))

    V = interpolated_V

    data[ind] =df_cd
    data_V[ind] = V
    
# Stack DataFrames and compute the mean
combined_df_mean = pd.concat(data.values()).groupby(level=0).mean()
combined_df_std = pd.concat(data.values()).groupby(level=0).std()

#predicted
data_p = {}
models = ['TPM', 'GM', 'WM_F0p5', 'SWM', 'UGM', 'UGM18']
for model in models:
    temp = {}
    for i, ind in enumerate(patientlist):
        df = pd.read_excel(f'predicted_cd/{model}_{ind}.xlsx', index_col = 0)
        temp[ind] = df
    df_mean = pd.concat(temp.values()).groupby(level=0).mean()
    df_std = pd.concat(temp.values()).groupby(level=0).std()
    data_p[model] = {'mean': df_mean, 'sd':df_std}
#%%
fig, axes = plt.subplots(2, 2, figsize=(8, 8), sharex=True)
axes = axes.flatten()

# Collect handles and labels only once to avoid repetition in the legend
handles, labels = [], []
model_names = ['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM18']
for i, solute in enumerate(df_cd.columns):
    # Plot the data with error bars
    h_data = axes[i].errorbar(combined_df_mean.index, combined_df_mean[solute], yerr=combined_df_std[solute], fmt='o', color='k', label='Data', capsize=5)
    
    # Collect handles and labels from the first plot
    if i == 0:
        handles.append(h_data)
        labels.append('Data')
    
    # Plot model predictions
    for j, model in enumerate(models):
        h_model, = axes[i].plot(np.arange(0, 241), data_p[model]['mean'][solute], ls='-', lw=2, label=model_names[j])
        
        # Collect handles and labels from the first plot
        if i == 0:
            handles.append(h_model)
            labels.append(model_names[j])
    
    # Set ticks and labels
    axes[i].set_xticks([60, 120, 180, 240])
    axes[i].set_xticklabels([1, 2, 3, 4], fontsize=12)
    axes[i].set_title(solute, fontsize=14)
    axes[i].tick_params(axis='both', labelsize=14)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['right'].set_visible(False)

# Set common x-axis labels
axes[2].set_xlabel('Time (hr)', fontsize=14)
axes[3].set_xlabel('Time (hr)', fontsize=14)

# Add a common y-axis label
fig.supylabel('Dialysate conc. (mmol/L)', fontsize=14)

# Add the legend below the subplots with 7 columns
fig.legend(handles=handles[:7], labels=labels[:7], loc='lower center', ncol=7, frameon=False)

# Adjust layout to make space for the legend
plt.tight_layout(rect=[0, 0.1, 1, 1])  # Adjust the bottom to fit the legend

# Save the figure
plt.savefig('Figure3.svg', dpi=600)

#%%

# fig, axes = plt.subplots(5, 4, figsize = (12, 8), sharex = True)
# axes = axes.flatten()

# for i,pfile in enumerate(patientlist):
#     predicted_cd = pd.read_excel(f'predicted_cd/TPM_{pfile}.xlsx', index_col = 0)
#     axes[i].plot(data[pfile].index, data[pfile]['Sodium'], marker = 'o', ls = '')
#     axes[i].plot(predicted_cd.index, predicted_cd['Sodium'])
#     axes[i].set_title(pfile)
# plt.savefig('final-plots/sodium-dialysate_concentration_human.png', dpi = 600)

# fig, axes = plt.subplots(5, 4, figsize = (12, 8), sharex = True)
# axes = axes.flatten()

# for i,pfile in enumerate(patientlist):
#     predicted_cd = pd.read_excel(f'predicted_cd/TPM_{pfile}.xlsx', index_col = 0)
#     axes[i].plot(data[pfile].index, data[pfile]['Glucose'], marker = 'o', ls = '')
#     axes[i].plot(predicted_cd.index, predicted_cd['Glucose'])
#     axes[i].set_title(pfile)
# plt.savefig('final-plots/glucose-dialysate_concentration_human.png', dpi = 600)

# fig, axes = plt.subplots(5, 4, figsize = (12, 8), sharex = True)
# axes = axes.flatten()

# for i,pfile in enumerate(patientlist):
#     predicted_cd = pd.read_excel(f'predicted_V/TPM_{pfile}.xlsx', index_col = 0)
#     axes[i].plot(np.arange(241), data_V[pfile], marker = 'o', ls = '')
#     axes[i].plot(np.arange(241), predicted_cd)
#     axes[i].set_title(pfile)
# plt.savefig('final-plots/glucose-dialysate_concentration_human.png', dpi = 600)

