# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 14:11:07 2023

@author: P70073624
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Assuming df_p1, df_p2, df_p3 are your dataframes
patient = 'P7'
df_p1 = pd.read_excel('predicted_cd/TPM_'+patient+'.csv.xlsx', index_col = 0)
df_p2 = pd.read_excel('predicted_cd/UGM_'+patient+'.csv.xlsx', index_col = 0)
df_p3 = pd.read_excel('predicted_cd/UGM18_'+patient+'.csv.xlsx', index_col = 0)
df_p4 = pd.read_excel('predicted_cd/GM_'+patient+'.csv.xlsx', index_col = 0)
df_p5 = pd.read_excel('predicted_cd/WM_F0p5_'+patient+'.csv.xlsx', index_col = 0)
df_p6 = pd.read_excel('predicted_cd/SWM_'+patient+'.csv.xlsx', index_col = 0)

labels = [ 'TPM','GM', 'WM', 'SWM', 'UGM', 'UGM-18']
#insert others here
dataframes = [df_p1,  df_p4, df_p5, df_p6, df_p2, df_p3]

# Plotting
solute_columns = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
num_rows = 2
num_cols = 4


#get experimental data
pfile = 'patient_files/pig1/session2/'+patient+'.csv'
df = pd.read_csv(pfile,skiprows = range(0,16), delimiter = "," \
                 , encoding= 'unicode_escape')
    
'''dialysate solute concentration'''
df_cd = pd.DataFrame(columns = solute_columns,dtype = float)
df_cd = df[solute_columns].iloc[10:18].copy()  #dialysate concentration
for column in df_cd:
    df_cd[column] = df_cd[column].astype('float64')  
df_cd.loc[:, "Creatinine"]  *= 0.001   
index = pd.Series([0,10,20,30,60,120,180, 240])
df_cd = df_cd.set_index([index])
#using .values here to copy only the values, otherwise it tries to match the indices of df and df_cd and it doesnt work
df_cd = df_cd.interpolate(method = 'index', limit_direction = "both")

'''dialysate volume'''
V = pd.read_csv(pfile,skiprows = range(0,45), delimiter = ",", \
                    encoding= 'unicode_escape')[["IP volume T=0 (mL)","IP volume T=240 (mL)"]].iloc[0] # IPV measured from haemoglobin
# print(df_V)
df_V = pd.read_csv(pfile,skiprows = range(0,60), delimiter = ",", \
                    encoding= 'unicode_escape')[["Time","IPV"]].iloc[2:10]
    

df_V['IPV'] = df_V['IPV'].replace(["#REF!",'#DIV/0!', '#VALUE!'], np.nan).interpolate()
df_V =df_V.astype('float64')
df_V.set_index('Time', inplace = True)
if np.isnan(df_V.iloc[0].iloc[0]):
    df_V.loc[0] = float(V.iloc[0])
    
if np.isnan(df_V.iloc[-1].iloc[0]):
    df_V.iloc[-1] = float(V.iloc[1])
    
df_V['IPV'] = df_V['IPV'].interpolate()


#Linear interpolation to find values of V at all times
f_V = interp1d(df_V.index, df_V, axis = 0)
interpolated_V = f_V(range(0,241))
V = interpolated_V.flatten()

# Set up the figure with a custom layout
fig = plt.figure(figsize=(15, 12))

# First row: 3 columns
ax1 = plt.subplot2grid((3, 3), (0, 0))
ax2 = plt.subplot2grid((3, 3), (0, 1))
ax3 = plt.subplot2grid((3, 3), (0, 2))

# Second row: 3 columns
ax4 = plt.subplot2grid((3, 3), (1, 0))
ax5 = plt.subplot2grid((3, 3), (1, 1))
ax6 = plt.subplot2grid((3, 3), (1, 2))

# Third row: Left subplot spans two columns, right subplot occupies one column
ax7 = plt.subplot2grid((3, 3), (2, 0), colspan=2)
ax8 = plt.subplot2grid((3, 3), (2, 2))

axes = [ax1, ax2, ax3, ax4, ax5, ax6]
for i, solute in enumerate(solute_columns):
    

    axes[i].plot(df_cd.index, df_cd[solute], label = 'data',marker = 'o', c= 'k', ls = '')

    x_pos = np.arange(241)
    for j, df in enumerate(dataframes):
        print(df)
        axes[i].plot(x_pos, df[solute], label = labels[j] )
    axes[i].set_title(solute, fontsize = 18)
    axes[i].set_xticks([60, 120, 180, 240])
    axes[i].set_xticklabels([1, 2, 3, 4])
    axes[i].tick_params('both', labelsize = 18)
axes[0].set_ylabel('Concentration (mmol/L)', fontsize = 18)
axes[3].set_ylabel('Concentration (mmol/L)', fontsize = 18)
axes[4].legend(fontsize = 12, frameon = False)
for ax in [ax4, ax5, ax6]:
    ax.set_xlabel('Time (hrs)', fontsize = 18)


#plot volume
predicted_V = pd.read_excel('predicted_V/TPM_'+patient+'.csv.xlsx', index_col = 0)
ax8.scatter(x_pos, V, c= 'k', label = 'data')
ax8.plot(x_pos, predicted_V, c = 'red', label =  'TPM predicted')
ax8.set_ylabel('V (ml)', fontsize = 18)
ax8.set_xlabel('Time(hrs)', fontsize = 18)
ax8.set_xticks([60, 120, 180, 240])
ax8.set_xticklabels([1, 2, 3, 4])
ax8.tick_params('both', labelsize = 18)
ax8.legend(fontsize = 10, frameon = False)
  
bar_width = 0.1
labels = [ 'TPM','GM', 'WM_F0p5', 'SWM', 'UGM', 'UGM18']
#plot MTAC
for i, model in enumerate(labels):
    mtac = pd.read_excel(f'fittedPS/fittedPS_{model}.xlsx', index_col = 0)  
    MTAC = mtac.loc[f'{patient}.csv'].iloc[:6]
    ax7.bar(np.arange(6)+i*bar_width, MTAC, bar_width, edgecolor = 'k', label = model)
ax7.set_ylabel('MTAC (ml/min)', fontsize = 18)
ax7.legend(frameon = False, ncols = 2)
ax7.set_xticks(np.arange(6)+0.25)
ax7.set_xticklabels(solute_columns)
ax7.tick_params('both', labelsize = 18)

label_positions = [-0.08, 1.1]  # Adjust as necessary
for i, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]):
    ax.text(label_positions[0], label_positions[1], chr(97 + i), transform=ax.transAxes, 
            fontsize=16, fontweight='bold', va='top', ha='right')
    
plt.tight_layout()
plt.show()
plt.savefig('final-plots/predicted_cd_'+patient+'.png', dpi = 600)
