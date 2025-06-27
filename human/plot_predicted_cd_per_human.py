# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:36:53 2024

@author: P70073624
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:55:41 2024

@author: P70073624
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

patientlist = ['PDPC_059', 'PDPC_236','PDPC_268','PDPC_541','PDPC_572','PDPC_800','PDPC_835','PDPC_890','PDPC_946']
patient_id = '800'
num_rows = len(patient_id)
# Set up the figure with a custom layout
fig = plt.figure(figsize=(18, 10))

# First row: 4 columns
ax1 = plt.subplot2grid((2, 4), (0, 0))
ax2 = plt.subplot2grid((2, 4), (0, 1))
ax3 = plt.subplot2grid((2, 4), (0, 2))
ax4 = plt.subplot2grid((2, 4), (0, 3))

# Second row: 3 columns and 1 column
ax5 = plt.subplot2grid((2, 4), (1, 0), colspan=3)
ax6 = plt.subplot2grid((2, 4), (1, 3))

axes = [ax1, ax2, ax3, ax4]


# Assuming df_p1, df_p2, df_p3 are your dataframes
df_p1 = pd.read_excel('predicted_cd/TPM_PDPC_'+patient_id+'.xlsx', index_col = 0)
df_p2 = pd.read_excel('predicted_cd/UGM_PDPC_'+patient_id+'.xlsx', index_col = 0)
df_p3 = pd.read_excel('predicted_cd/UGM18_PDPC_'+patient_id+'.xlsx', index_col = 0)
df_p4 = pd.read_excel('predicted_cd/GM_PDPC_'+patient_id+'.xlsx', index_col = 0)
df_p5 = pd.read_excel('predicted_cd/WM_F0p5_PDPC_'+patient_id+'.xlsx', index_col = 0)
df_p6 = pd.read_excel('predicted_cd/SWM_PDPC_'+patient_id+'.xlsx', index_col = 0)

labels = ['TPM',  'GM', 'WM', 'SWM', 'UGM', 'UGM-18']
#insert others here
dataframes = [df_p1,  df_p4, df_p5, df_p6, df_p2, df_p3]

# Plotting
solute_columns = df_p1.columns



#get experimental data

datafile = pd.read_excel('PD_in_silico_csv_export_/PD_in_silico_export_20240723.xlsx')
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)
ind = 'PDPC_'+patient_id

data = {}


# dialysate concentrations
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
tp_ond = datafile.loc[ind, "TP_OND"] #total protein overnight dwell in g/L
tp_ff = datafile.loc[ind, "TP_FF"]

if np.isnan(tp_ond) or tp_ond < 1:
    tp_ond = datafile["TP_OND"].mean()

if np.isnan(tp_ff):
    tp_ff = datafile["TP_FF"].mean()
    
# RV = Vfill*tp_ff/(tp_ond-tp_ff)
RV0 = datafile.loc[ind, "RV0"]
RV240 = datafile.loc[ind, "RV240"]

ur_0 = [datafile.loc[ind, 'Ur_night']*RV0/(Vfill+RV0) if datafile.loc[ind, 'Ur_night'] > 0 else 0]
cr_0 = [datafile.loc[ind, 'Cr_night']*RV0/(Vfill+RV0)*0.001 if datafile.loc[ind, 'Cr_night'] > 0 else 0]
glu_0 = [214 if datafile.loc[ind, 'Gluc_PET'] == 3.86 else  126 if datafile.loc[ind, 'Gluc_PET'] == 2.27 else  75 ]
na_0 = [132]
df_cd.loc[0] = np.concatenate((ur_0, cr_0, glu_0, na_0))


# Volume from datafiles
df_V = np.array(datafile.loc[ind,['Vin_t0', 'Vdrain_t240']]) #mL
df_V[0] += RV0 #assuming a residual volume of 200 ml
df_V[1] += RV240
# df_V += RV #assuming a residual volume of 200 ml
#Linear interpolation to find values of V at all times
f_V = interp1d([0,240], df_V)
interpolated_V = f_V(range(0,241))

V = interpolated_V
    



for i, solute in enumerate(solute_columns):
    

    axes[i].plot(df_cd.index, df_cd[solute], label = 'data',marker = 'o', c= 'k', ls = '')

    x_pos = np.arange(241)
    for j, df in enumerate(dataframes):
        axes[i].plot(x_pos, df[solute], label = labels[j], lw =3 )
        axes[i].set_title(solute, fontsize = 18)
        axes[i].tick_params(axis='both', labelsize=18)
        axes[i].set_xlabel('Time (hrs)', fontsize = 18)
        axes[i].set_xticks([60, 120, 180, 240])
        axes[i].set_xticklabels([1, 2, 3, 4])
        # axes[0].text(0.1, 0.9, ind[0], transform = axes[0].transAxes)
axes[2].legend( frameon = False)     
axes[0].set_ylabel('Concentration (mmol/L)', fontsize = 18)
   
#plot volume
predicted_V = pd.read_excel('predicted_V/TPM_PDPC_'+patient_id+'.xlsx', index_col = 0)
ax6.scatter(x_pos, V, c= 'k', label = 'data')
ax6.plot(x_pos, predicted_V,  label =  'TPM predicted')
ax6.set_ylabel('V (ml)', fontsize = 18)
ax6.set_xlabel('Time(hrs)', fontsize = 18)
ax6.set_xticks([60, 120, 180, 240])
ax6.set_xticklabels([1, 2, 3, 4])
ax6.legend( frameon = False)
ax6.tick_params('both', labelsize = 18)

bar_width = 0.1
labels = [ 'TPM','GM', 'WM_F0p5', 'SWM', 'UGM', 'UGM18']
#plot MTAC
for i, model in enumerate(labels):
    mtac = pd.read_excel(f'fittedPS/fittedPS_{model}.xlsx', index_col = 0)  
    MTAC = mtac.loc[f'PDPC_{patient_id}'].iloc[:4]
    ax5.bar(np.arange(4)+i*bar_width, MTAC, bar_width, edgecolor = 'k', label = model)
ax5.set_ylabel('MTAC (ml/min)', fontsize = 18)
ax5.legend(frameon = False, ncols = 2)
ax5.set_xticks(np.arange(4)+0.15)
ax5.set_xticklabels(solute_columns, fontsize = 18)
ax5.tick_params('y', labelsize = 18)

label_positions = [-0.08, 1.1]  # Adjust as necessary
for i, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6]):
    ax.text(label_positions[0], label_positions[1], chr(97 + i), transform=ax.transAxes, 
            fontsize=16, fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.show()
plt.savefig(f'final-plots/predicted_cd_PDPC_{patient_id}.png', dpi = 600)