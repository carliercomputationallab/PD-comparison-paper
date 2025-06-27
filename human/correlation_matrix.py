# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 14:45:58 2024
Draw correlation between demographics and fitted parameters MTACs and 
hydraulic conductivity
@author: P70073624
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# read patient data
datafile = pd.read_csv('PD_in_silico_csv_export_/PD_in_silico_export_20240723.csv', sep = ";")
datafile.set_index('Participant Id', inplace=True)
patientlist = datafile.index


# append MTAC and LPS -  fitted parameters to the data frame
MTAC = pd.read_excel('fittedPS/fittedPS_TPM.xlsx', index_col=0)
MTAC = MTAC[['MTAC_urea', 'MTAC_crea', 'MTAC_glu', 'MTAC_sodium', 'L', 'total_ultrafiltration (ml)']]
MTAC['Ultrafiltration'] = MTAC['total_ultrafiltration (ml)'] + MTAC['L']*240
result = pd.concat([datafile, MTAC], axis=1)
result = result.drop(result.columns[[0, 1, 2]], axis=1)
result.drop(['Vin_night', 'Vdrain_night','Cr_night', 'Ur_night','Creatinine_t20', 'Creatinine_t120',
'Creatinine_t240', 'Gluc_t20', 'Gluc_t120', 'Gluc_t240', 'Ur_t20',
'Ur_t120', 'Ur_t240', 'Na20', 'Na120', 'Na240',
'TP_OND', 'TP_FF', 'TP_t240', 'TP_NFF', 'total_ultrafiltration (ml)'], axis = 1, inplace = True)
result = result.rename(columns = {'Vin_t0': 'Fill volume', 'RV_t0': 'Residual volume', 'Cr_t2_B': 'Cp, creatinine',
                                  'Ur_t2_B':'Cp, urea', 'Gluc_t2_B':'Cp, glucose', 'Na_t2_B': 'Cp, Sodium',
                                  "Residual_Diuresis":"Residual Diuresis", "Gluc_PET": "Glucose PET",
                                  "Vdrain_t240": "Drain volume", "MTAC_urea": "MTAC, urea", "MTAC_crea":"MTAC, creatinine",
                                  "MTAC_glu": "MTAC, glucose", "MTAC_sodium":"MTAC, sodium", 'L': 'Lymphatic flow rate', 'Sex': 'Sex (F)'})
corr = result.corr()

# Create a mask for values greater than 0.5 or less than -0.5
mask = (corr.abs() > 0.5) | (corr.abs() < -0.5)

# Set up the figure and axis
fig, ax = plt.subplots(figsize=(15, 25))


# Plot the matrix with only values satisfying the condition
cax = ax.matshow(np.tril(corr[mask]), cmap='coolwarm')

# Add x ticks and y ticks with column names
ax.set_xticks(range(len(result.columns)))
ax.set_yticks(range(len(result.columns)))

# Prevent x-tick labels from overlapping
ax.tick_params(axis='x', rotation=90, labelsize=18)

# Set half of the y-tick labels to a different font color and bold
for i, label in enumerate(ax.get_yticklabels()):
    if i < 14:
        label.set_color('red')
        # label.set_weight('bold')
        
# Set half of the y-tick labels to a different font color and bold
for i, label in enumerate(ax.get_xticklabels()):
    if i < 14:
        label.set_color('red')
        # label.set_weight('bold')

ax.tick_params(axis='y', labelsize=18)  # Adjust y-tick label size if needed

# Set x-tick labels (possibly adjusting fontsize)
ax.set_xticklabels(result.columns, rotation=90, fontsize=18)

# Set y-tick labels (possibly adjusting fontsize)
ax.set_yticklabels(result.columns, fontsize=18)

# Display the colorbar
plt.colorbar(cax)

# Adjust layout to prevent overlap of x-tick labels
plt.tight_layout()
# Manually draw grid lines centered on cells by offsetting by 0.5
for i in range(len(result.columns) + 1):
    ax.axhline(i - 0.5, color='grey', linestyle='-', linewidth=0.5)
    ax.axvline(i - 0.5, color='grey', linestyle='-', linewidth=0.5)

# Ensure the grid lines are behind the matshow plot
ax.set_axisbelow(True)


# Show the plot
plt.show()

plt.savefig('demographics-correlation_matrix.png', dpi = 600)


