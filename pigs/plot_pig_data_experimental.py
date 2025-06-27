# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:49:27 2024
plot patient data
@author: P70073624
"""



from values import input_values
import os
import matplotlib.pyplot as plt
from fnmatch import fnmatch

"Get all MTACs for the pig in question"
root = 'patient_files/'
pattern = "*.csv"
patientlist = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            patientlist.append(os.path.join(path, name))
# fig, ax = plt.subplots(3,2, figsize = (16,18))
# for pfile in patientlist:
#     predicted_cd, Cp, V, df_cd, Vr, V_fill = input_values(pfile)
#     V = V*1000
#     t = 240
#     #urea
#     df_cd['Urea'].plot( ax = ax[0,0], label = 'data', style = '.')
#     ax[0,0].set_title("Urea")
    
#     #creatinine
#     df_cd['Creatinine'].plot( ax = ax[0,1], style = '.')
#     ax[0,1].set_title("Creatinine")
    
#     #Sodium
#     df_cd['Sodium'].plot( ax = ax[1,0],  style = '.')
#     ax[1,0].set_title("Sodium")
    
#     #Phosphate
#     df_cd['Phosphate'].plot( ax = ax[1,1], style = '.')
#     ax[1,1].set_title("Phosphate")
    
#     #Glucose
#     df_cd['Glucose'].plot( ax = ax[2,0], style = '.')
#     ax[2,0].set_title("Glucose")
    
#     #Potassium
#     df_cd['Potassium'].plot( ax = ax[2,1], style = '.', label = f'{pfile[31:]}')
#     ax[2,1].set_title("Potassium")
    
# fig.supxlabel("time, min")
# fig.supylabel("Dialysate concentration, mmol")
# handles, labels = ax[2, 1].get_legend_handles_labels()
# fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)

# plt.subplots_adjust(top=0.88,
#                     bottom=0.11,
#                     left=0.09,
#                     right=0.9,
#                     hspace=0.295,
#                     wspace=0.215)


fig, ax = plt.subplots(6,5, figsize = (16,18), sharex = True)
axes = ax.flatten()
solute = 'Sodium'
for i,pfile in enumerate(patientlist):
    predicted_cd, Cp, V, df_cd, Vr, V_fill = input_values(pfile)
    V = V*1000
    t = 240
    
    ax2 = axes[i].twinx()
    #Sodium
    df_cd[solute].plot( ax = axes[i],  style = '.-', label = 'dial')
    axes[i].set_title( f'{pfile[28:]}')
    ax2.plot(range(t+1),Cp[:,2],  ls = '-', color = 'red', label = 'blood')  
    # Format the secondary y-axis
    ax2.spines['right'].set_color('red')
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red')
axes[i].legend(loc = 'upper left')    
ax2.legend()   
fig.supxlabel("time, min")
fig.supylabel(f"Dialysate {solute} concentration, mmol")


plt.subplots_adjust(top=0.975,
                    bottom=0.11,
                    left=0.08,
                    right=0.965,
                    hspace=0.48,
                    wspace=0.5)
