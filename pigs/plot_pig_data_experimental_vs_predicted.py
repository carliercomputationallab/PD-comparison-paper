# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 14:13:30 2024

@author: P70073624
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from values import input_values
from fnmatch import fnmatch
import os

"Get all MTACs for the pig in question"
root = 'patient_files/'
pattern = "*.csv"
patientlist = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            patientlist.append(os.path.join(path, name))

terms_to_remove = ['P1.csv', 'P14.csv', 'P15.csv', 'P18.csv', 'P19.csv', 'P22.csv', 'P23.csv', 'P28.csv', 'P29.csv']
data = {}
cp = {}

for i,pfile in enumerate(patientlist):
    ind = pfile[28:]
    if ind not in terms_to_remove:
        predicted_cd, Cp, V, df_cd, Vr, V_fill,_ = input_values(pfile)
        data[pfile[28:-4]] = df_cd
        cp[pfile[28:-4]] = Cp
    
# Stack DataFrames and compute the mean
combined_df_mean = pd.concat(data.values()).groupby(level=0).mean()
combined_df_std = pd.concat(data.values()).groupby(level=0).std()

#predicted
data_p = {}
models = ['TPM', 'GM', 'WM_F0', 'SWM', 'UGM', 'UGM18']

for model in models:
    temp = {}
    for i, pfile in enumerate(patientlist):
        ind = pfile[28:]
        if ind not in terms_to_remove:
            df = pd.read_excel(f'predicted_cd/{model}_{ind}.xlsx', index_col = 0)
            temp[ind] = df
    df_mean = pd.concat(temp.values()).groupby(level=0).mean()
    df_std = pd.concat(temp.values()).groupby(level=0).std()
    data_p[model] = {'mean': df_mean, 'sd':df_std}
#%% 
fig, axes = plt.subplots(3, 2, figsize=(8, 9), sharex='col')
axes = axes.flatten()
model_names = ['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM18']
# Initialize an empty list to collect handles and labels
handles, labels = [], []

for i, solute in enumerate(df_cd.columns):
    # Plot the data with error bars
    h_data = axes[i].errorbar(combined_df_mean.index, combined_df_mean[solute], yerr=combined_df_std[solute], fmt='o', color='k', label='Data', capsize=5)
    
    # Collect the handle and label only for the first plot
    if i == 0:
        handles.append(h_data)
        labels.append('Data')
    
    # Plot the model predictions
    for j, model in enumerate(models):
        h_model, = axes[i].plot(np.arange(241), data_p[model]['mean'][solute], ls='-', lw=2, label=model_names[j])
        
        # Collect handles and labels only for the first plot
        if i == 0:
            handles.append(h_model)
            labels.append(model_names[j])
    
    # Set ticks and labels
    axes[i].set_xticks([60, 120, 180, 240])
    axes[i].set_xticklabels([1, 2, 3, 4])
    axes[i].set_title(solute, fontsize=14)
    axes[i].tick_params(axis='x', labelsize=12)
    axes[i].tick_params(axis='y', labelsize=12)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['right'].set_visible(False)

# Set common x-axis labels
axes[4].set_xlabel('Time (hr)', fontsize=14)
axes[5].set_xlabel('Time (hr)', fontsize=14)

# Add a common y-axis label
fig.supylabel('Dialysate conc. (mmol/L)', fontsize=14)

# Add the legend using only the first set of handles and labels, with 7 columns, centered
fig.legend(handles=handles[:7], labels=labels[:7], loc='lower center', ncol=7, frameon=False)

# Adjust layout to make space for the legend
# plt.tight_layout(rect=[0, 0.05, 1, 1])

# Save the figure
plt.savefig('Figure2.svg', dpi=600)

#%%
fig, axes = plt.subplots(5, 6, figsize = (12, 8), sharex = True)
axes = axes.flatten()
for i,pfile in enumerate(patientlist):
    ind = pfile[28:]
    if ind not in terms_to_remove:
        predicted_cd, Cp, V, df_cd, Vr, V_fill,_ = input_values(pfile)
        predicted_cd = pd.read_excel(f'predicted_cd/TPM_{ind}.xlsx', index_col = 0)
        axes[i].plot(df_cd.index, df_cd['Sodium'], marker = 'o', ls = '')
        axes[i].plot(predicted_cd.index, predicted_cd['Sodium'])
    axes[i].set_title(ind)
plt.savefig('final-plots/sodium-dialysate_concentration_pigs.png', dpi = 600)
fig, axes = plt.subplots(5, 6, figsize = (12, 8), sharex = True)
axes = axes.flatten()
for i,pfile in enumerate(patientlist):
    ind = pfile[28:]
    if ind not in terms_to_remove:
        predicted_cd, Cp, V, df_cd, Vr, V_fill,_ = input_values(pfile)
        predicted_V = pd.read_excel(f'predicted_V/TPM_{ind}.xlsx', index_col = 0)
        axes[i].plot(np.arange(241), V, marker = 'o', ls = '')
        axes[i].plot(np.arange(241), predicted_V.iloc[:,0])
    axes[i].set_title(ind)
plt.savefig('final-plots/predicted_volume_pigs.png', dpi = 600)