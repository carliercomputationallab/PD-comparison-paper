# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 21:10:17 2024

@author: P70073624
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from fnmatch import fnmatch

# Data (as manually entered from your example)
df_pigs = pd.read_excel('pigs/fittedPS/fittedPS_TPM.xlsx', index_col = 0)
df_human  = pd.read_excel('human/fittedPS/fittedPS_TPM.xlsx', index_col = 0)
terms_to_remove = ['P1.csv', 'P14.csv', 'P15.csv', 'P18.csv', 'P19.csv', 'P22.csv', 'P23.csv', 'P28.csv', 'P29.csv'] #extraneal dwells we are focusing on physioneal
df_pigs.drop(labels=terms_to_remove, inplace = True)

# Optional: Reorder df_human columns to match df_pigs for consistency
df_human = df_human[['total_ultrafiltration (ml)', 'L', 'V_err']]
df_pigs = df_pigs[['Total UF', 'L', 'V_err']]

df_human.columns = ['Total UF', 'L', 'V_err']

df_pigs['UF'] = df_pigs['Total UF'] 
df_human['UF'] = df_human['Total UF'] 
# get pigs data
root = 'pigs/patient_files/'
pattern = "*.csv"
patientlist = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            patientlist.append(os.path.join(path, name))


patientlist = [pfile for pfile in patientlist 
                        if not any(term in pfile for term in terms_to_remove)]
glu0 = []    
glu_sol = []       
for pfile in patientlist:
    print(pfile)
    df = pd.read_csv(pfile,skiprows = range(0,16), delimiter = "," \
                     , encoding= 'unicode_escape')
    solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
    
    '''dialysate solute concentration'''
    df_cd = pd.DataFrame(columns = solutes,dtype = float)
    df_cd = df[solutes].iloc[10:18].copy()  #dialysate concentration
    for column in df_cd:
        df_cd[column] = df_cd[column].astype('float64')  
    df_cd.loc[:, "Creatinine"]  *= 0.001   
    index = pd.Series([0,10,20,30,60,120,180, 240])
    df_cd = df_cd.set_index([index])
    #using .values here to copy only the values, otherwise it tries to match the indices of df and df_cd and it doesnt work
    df_cd = df_cd.interpolate(method = 'index', limit_direction = "both")   
    cd = df_cd['Glucose'].loc[0]
    
    glu0.append(cd)
    if cd < 75:
        sol = 1.36
    elif cd> 214:
        sol = 3.86
    else:
        sol = 2.27
        
    glu_sol.append(sol)

df_pigs['glucose %'] = glu_sol


#get human data
datafile = pd.read_csv('human/PD_in_silico_csv_export_/PD_in_silico_export_20240723.csv', sep = ",")
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)
patientlist = ['PDPC_059', 'PDPC_236','PDPC_268','PDPC_541','PDPC_572','PDPC_800','PDPC_835','PDPC_890','PDPC_946']
glu0 = []
for ind in patientlist:
    glu0.append(datafile.loc[ind, 'Gluc_PET'])
df_human = df_human[df_human.index.isin(patientlist)]    
df_human['glucose %'] = glu0


df_pigs['Group'] = 'pigs'
df_human['Group'] = 'humans'
# Combine the two DataFrames into one for easier plotting
combined_df = pd.concat([df_pigs, df_human], ignore_index=True)
# Ensure 'glucose %' 3.86 is represented in 'pigs'
for glucose in combined_df['glucose %'].unique():
    if not ((combined_df['glucose %'] == glucose) & (combined_df['Group'] == 'pigs')).any():
        # Adding a row with NaN for 'Total UF' to represent the group
        # Create a DataFrame with the new row
        new_row = pd.DataFrame({'Total UF': [np.nan], 'glucose %': [glucose], 'Group': ['pigs']})
        # Concatenate with the existing DataFrame
        combined_df = pd.concat([combined_df, new_row], ignore_index=True)
combined_df.sort_values(['Group', 'glucose %'], ascending = False, inplace = True)
# Define conditions for each group number
conditions = [
    (combined_df['glucose %'] == 1.36) & (combined_df['Group'] == 'pigs'),
    (combined_df['glucose %'] == 1.36) & (combined_df['Group'] == 'humans'),
    (combined_df['glucose %'] == 2.27) & (combined_df['Group'] == 'pigs'),
    (combined_df['glucose %'] == 2.27) & (combined_df['Group'] == 'humans'),
    (combined_df['glucose %'] == 3.86) & (combined_df['Group'] == 'pigs'),
    (combined_df['glucose %'] == 3.86) & (combined_df['Group'] == 'humans')
]

# Define corresponding values for each condition
group_numbers = [-0.2, 0.2, 0.8, 1.2, 1.8, 2.2]

# Use np.select to assign group numbers based on conditions
combined_df['group number'] = np.select(conditions, group_numbers, default=np.nan)



# Literature values for L and V_err
Literature_L = 0.3
Literature_V_err = np.nan  # Replace with the actual value if available

# Set up the figure with a custom layout
fig = plt.figure(figsize=(10, 8))

# Left side subplot (Code 1)
ax_left = plt.subplot2grid((2, 2), (0, 0), rowspan=2)

# Right side subplots (Code 2)
ax_right_top = plt.subplot2grid((2, 2), (0, 1))
ax_right_bottom = plt.subplot2grid((2, 2), (1, 1))

# Plotting for Code 1 data (Total UF)
sns.violinplot(x='glucose %', y='UF', data=combined_df, ax=ax_left, hue='Group', palette=['#CA1F7B', '#0070BB'], legend=True)
sns.scatterplot(x='group number', y='UF', data=combined_df, ax=ax_left, hue='Group', palette=['#FF91AF', '#89CFF0'], s=100, edgecolor='w', linewidth=1, legend=False)
ax_left.legend(frameon=False, loc='lower right')
ax_left.set_ylabel('Total UF [ml]', fontsize=14)
ax_left.spines['top'].set_visible(False)
ax_left.spines['right'].set_visible(False)

# Adding label "a" to the left subplot
ax_left.text(-0.1, 1.05, 'a', transform=ax_left.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

# Plotting for variable L (Code 2) on the top right
sns.violinplot(x='Group', y='L', data=combined_df, hue='Group', inner='quartile', ax=ax_right_top, palette=['#CA1F7B', '#0070BB'])
ax_right_top.set_ylabel('L [ml/min]', fontsize=14)
ax_right_top.set_xticks([-0.2, 1.2])
ax_right_top.set_xticklabels('')
ax_right_top.set_xlabel('')
ax_right_top.legend().set_visible(False)
ax_right_top.tick_params(axis='y', labelsize=14)
ax_right_top.spines['top'].set_visible(False)
ax_right_top.spines['right'].set_visible(False)

# Adding label "b" to the top right subplot
ax_right_top.text(-0.1, 1.05, 'b', transform=ax_right_top.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')


# Plotting for variable V_err (Code 2) on the bottom right
sns.violinplot(x='Group', y='V_err', data=combined_df, hue='Group', inner='quartile', ax=ax_right_bottom, palette=['#CA1F7B', '#0070BB'])
ax_right_bottom.set_ylabel('Volume error', fontsize=14)
ax_right_bottom.set_xticks([-0.2, 1.2])
ax_right_bottom.set_xticklabels(['pigs', 'humans'])
ax_right_bottom.set_xlabel('')
ax_right_bottom.legend().set_visible(False)
ax_right_bottom.tick_params(axis='y', labelsize=14)
ax_right_bottom.spines['top'].set_visible(False)
ax_right_bottom.spines['right'].set_visible(False)

# Adding label "c" to the bottom right subplot
ax_right_bottom.text(-0.1, 1.05, 'c', transform=ax_right_bottom.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

# Adjust layout and show
plt.tight_layout()
plt.show()

# Save the combined plot
plt.savefig('combined_uf_verr_L.png', dpi=600)