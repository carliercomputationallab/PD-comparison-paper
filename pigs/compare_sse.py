# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 15:59:34 2023

@author: P70073624
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Assuming df_p1, df_p2, df_p3 are your dataframes
df_p1 = pd.read_excel('fittedPS/sse_m7.xlsx', index_col = 0)
df_p2 = pd.read_excel('fittedPS/sse_unifiedGraff.xlsx', index_col = 0)
df_p3 = pd.read_excel('fittedPS/sse_unifiedGraff_18p.xlsx', index_col = 0)
df_p4 = pd.read_excel('fittedPS/sse_m8.xlsx', index_col = 0)
df_p5 = pd.read_excel('fittedPS/sse_Waniewski.xlsx', index_col = 0)
df_p6 = pd.read_excel('fittedPS/sse_m9.xlsx', index_col = 0)
labels = ['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM-18']
#insert others here
dataframes = [df_p1, df_p4, df_p5, df_p6, df_p2, df_p3]

# Plotting
solute_columns = df_p1.columns
num_rows = 2
num_cols = 3

fig, axs = plt.subplots(num_rows, num_cols, figsize=(9, 6), sharex= True)

for i, solute in enumerate(solute_columns):
    row_index = i // num_cols
    col_index = i % num_cols
    
    means = [df[solute].mean() for df in dataframes]
    stds = [df[solute].std() for df in dataframes]
    

    x_pos = np.arange(len(labels))

    axs[row_index, col_index].errorbar(x_pos, means, yerr=stds, color = 'k', marker = 'o', alpha=0.7, ecolor='black', capsize=5, ls = '')
    axs[row_index, col_index].set_title(solute)
    
    for j, df in enumerate(dataframes):
        axs[row_index, col_index].scatter([x_pos[j]] * len(df), df[solute], s=10, alpha=0.5, color='gray')


    

axs[1, 0].set_xticks(x_pos)
axs[1, 0].set_xticklabels(labels, rotation = 30)
axs[1, 1].set_xticks(x_pos)
axs[1, 1].set_xticklabels(labels, rotation = 30)
axs[1, 2].set_xticks(x_pos)
axs[1, 2].set_xticklabels(labels, rotation = 30)
fig.supylabel('SSE')
plt.tight_layout()
plt.show()

plt.savefig('final-plots/sse.png', dpi = 600)
