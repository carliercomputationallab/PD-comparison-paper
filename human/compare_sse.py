# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 15:03:15 2024
compare solute specific error for each model
@author: P70073624
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# store all sse
df_p1 = pd.read_excel('fittedPS/sse_m7.xlsx', index_col = 0)
df_p2 = pd.read_excel('fittedPS/sse_unifiedGraff.xlsx', index_col = 0)
df_p3 = pd.read_excel('fittedPS/sse_unifiedGraff_18p.xlsx', index_col = 0)
df_p4 = pd.read_excel('fittedPS/sse_m8.xlsx', index_col = 0)
df_p5 = pd.read_excel('fittedPS/sse_Waniewski.xlsx', index_col = 0)
df_p6 = pd.read_excel('fittedPS/sse_m9.xlsx', index_col = 0)
labels = ['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM-18']

dataframes = [df_p1, df_p4, df_p5, df_p6, df_p2, df_p3]

# Plotting
solute_columns = df_p1.columns
num_rows = 1
num_cols = 3

fig, axs = plt.subplots(num_rows, num_cols, figsize=(9, 3), sharex= True)

for i, solute in enumerate(solute_columns):

    
    means = [df[solute].mean() for df in dataframes]
    stds = [df[solute].std() for df in dataframes]
    

    x_pos = np.arange(len(labels))

    axs[i].errorbar(x_pos, means, yerr=stds, color = 'k', marker = 'o', alpha=0.7, ecolor='black', capsize=5, ls = '')
    axs[i].set_title(solute)
    
    # for j, df in enumerate(dataframes):
    #     axs[i].scatter([x_pos[j]] * len(df), df[solute], s=10, alpha=0.5, color='gray')


    

axs[0].set_xticks(x_pos)
axs[0].set_xticklabels(labels, rotation = 30)
axs[1].set_xticks(x_pos)
axs[1].set_xticklabels(labels, rotation = 30)
axs[2].set_xticks(x_pos)
axs[2].set_xticklabels(labels, rotation = 30)
fig.supylabel('SSE')
plt.tight_layout()
plt.show()

plt.savefig('final-plots/sse.png', dpi = 600)
