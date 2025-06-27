import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Data (as manually entered from your example)
df_pigs = pd.read_excel('pigs/fittedPS/fittedPS_TPM.xlsx', index_col = 0)
df_human  = pd.read_excel('human/fittedPS/fittedPS_TPM.xlsx', index_col = 0)
cols = ['MTAC_urea', 'MTAC_crea', 'MTAC_sodium', 'MTAC_phosphate', 'MTAC_glu',
       'MTAC_potassium', 'L', 'V_err']
for column in cols:
    if column not in df_human.columns:
        df_human[column] = np.nan

# Optional: Reorder df_human columns to match df_pigs for consistency
df_human = df_human[cols]
df_pigs = df_pigs[cols]
df_human.columns = ['Urea', 'Creatinine', 'Sodium', 'Phosphate', 'Glucose',
       'Potassium', 'L', 'V_err']
df_pigs.columns = ['Urea', 'Creatinine', 'Sodium', 'Phosphate', 'Glucose',
       'Potassium', 'L', 'V_err']
df_pigs['Group'] = 'pigs'
df_human['Group'] = 'human'
# Combine the two DataFrames into one for easier plotting
combined_df = pd.concat([df_pigs, df_human], ignore_index=True)
#%%
# MTAC and LpS values from Oberg 2019 paper
Literature = [15.11, 5.682, 1.936, 4.171, 7.319, 24.352, 0.3, np.nan]

# Set up the figure and subplots
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(10, 8))  # One row, eight columns
axes = axes.flatten()

for i, col in enumerate(df_pigs.columns[:-1]):  # Ignore the 'Group' column
    sns.violinplot(x='Group', y=col, data=combined_df, hue='Group', inner='quartile', ax=axes[i], palette=['#CA1F7B', '#0070BB'])
    axes[i].plot(1, Literature[i], ls = '', marker = 'o', color = '#C41E3A',mec = 'yellow', ms = 8, label = 'TPM lit.')
    axes[i].set_ylabel(f'MTAC {col} [ml/min]', fontsize = 16)
    axes[i].set_xticklabels('')
    axes[i].set_xlabel('')
    axes[i].legend().set_visible(False)
    axes[i].tick_params(axis = 'y', labelsize =16)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['right'].set_visible(False)
axes[0].legend(loc = 'upper left', frameon = False)
axes[6].set_ylabel('L [ml/min]', fontsize = 16)
axes[7].set_ylabel('Volume error', size = 16)

# Adjust layout
plt.tight_layout()
plt.show()

plt.savefig('compare_fitted_parameters.png', dpi = 600)
