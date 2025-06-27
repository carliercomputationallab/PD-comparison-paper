# -s- coding: utf-8 -s-
"""
Created on Fri Apr  5 17:15:44 2024
compare_RMSE_pigs_human
@author: P70073624
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines

from PIL import Image
# Open an image from a computer 
def open_image_local(path_to_image):
    image = Image.open(path_to_image) # Open the image
    image_array = np.array(image) # Convert to a numpy array
    return image_array # Output


human_folder = 'human/fittedPS/'
# Read the DataFrame from a file (adjust the file path and format accordingly)
df_h1= pd.read_excel(human_folder+'fittedPS_TPM.xlsx')
df_h2 = pd.read_excel(human_folder+'fittedPS_GM.xlsx')
df_h3 = pd.read_excel(human_folder+'fittedPS_WM_F0p5.xlsx')
# df_W0 = pd.read_excel(human_folder+'fittedPS_WM_F0.xlsx')
# df_W1 = pd.read_excel(human_folder+'fittedPS_WM_F1.xlsx')
df_h4 = pd.read_excel(human_folder+'fittedPS_SWM.xlsx')
df_h5 = pd.read_excel(human_folder+'fittedPS_UGM.xlsx')
df_h6 = pd.read_excel(human_folder+'fittedPS_UGM18.xlsx')


df_list = [df_h1, df_h2, df_h3, df_h4, df_h5, df_h6]

###############################################
# PLOT OBJECTIVE FN AND COMP TIME
###############################################

# Initialize lists to store mean and std values
mean_obj_fn = []
std_obj_fn = []
mean_comp_time = []
std_comp_time = []

# Compute mean and std for each dataframe
for i, df in enumerate(df_list):
    if i == 0:  # Use C_err for the first model
        mean_obj_fn.append(df['C_err'].mean())
        std_obj_fn.append(df['C_err'].std())
    else:
        mean_obj_fn.append(df['obj_fn'].mean())
        std_obj_fn.append(df['obj_fn'].std())
    mean_comp_time.append(df['comp-time'].mean())
    std_comp_time.append(df['comp-time'].std())
print(mean_obj_fn, std_obj_fn)
# Create a figure
fig, ax = plt.subplots(1, 2, figsize=(8,4))
ax1 = ax[0]
ax2 = ax[1]
# Plot mean and std of 'obj_fn' on the left y-axis
x_values = np.arange(len(df_list))
ax1.errorbar(x_values, mean_obj_fn, yerr=std_obj_fn,  marker='o', color='#0070BB', capsize=5, ls = '')
# mean_W0 = df_W0['obj_fn'].mean()
# std_W0 = df_W0['obj_fn'].std()
# ax1.errorbar(1.8, mean_W0, yerr= std_W0, marker='d', color='#003262', capsize=5, ls = '' )
# mean_W1 = df_W1['obj_fn'].mean()
# std_W1 = df_W1['obj_fn'].std()
# ax1.errorbar(2.2, mean_W1, yerr= std_W1, marker='d', color='#7BAFD4', capsize=5, ls = '' )
ax1.set_ylabel(r'C$_{error}$', color='k')
# ax1.tick_params('y', colors='b')
ax1.set_xticks(x_values)
ax1.set_xticklabels(['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM-18'], rotation = 45)
ax1.text(0.05, 0.95, '(a)', transform = ax1.transAxes)
# Create a secondary y-axis for 'comp_time'
# ax2 = ax1.twinx()
ax2.errorbar(x_values,np.array(mean_comp_time) / 1000, yerr=np.array(std_comp_time) / 1000, marker='o', color='#0070BB', capsize=5, ls = '')
ax2.set_ylabel('Computational time, x 1000 s')
ax2.set_xticks(x_values)
ax2.set_xticklabels(['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM-18'], rotation = 45)
ax2.text(0.05, 0.95, '(b)', transform = ax2.transAxes)
# ax2.tick_params('y', colors='r')
# plt.subplots_adjust(top=0.94,
#                     bottom=0.23,
#                     left=0.125,
#                     right=0.9,
#                     hspace=0.2,
#                     wspace=0.2)

####################################################
#############              PIGS          ###########
####################################################
pig_folder = 'pigs/'
# Read the DataFrame from a file (adjust the file path and format accordingly)
df_p1 = pd.read_excel(pig_folder +'fittedPS/fittedPS_TPM.xlsx')
df_p2 = pd.read_excel(pig_folder +'fittedPS/fittedPS_GM.xlsx')
# Replace all missing values in the entire DataFrame with NaN
df_p2.replace(to_replace='missing_value', value=float('nan'), inplace=True)

df_p3 = pd.read_excel(pig_folder +'fittedPS/fittedPS_WM_F0p5.xlsx')
# df_W0 = pd.read_excel(pig_folder +'fittedPS/fittedPS_WM_F0.xlsx')
# df_W1 = pd.read_excel(pig_folder +'fittedPS/fittedPS_WM_F1.xlsx')
df_p4 = pd.read_excel(pig_folder +'fittedPS/fittedPS_SWM.xlsx')
df_p5 = pd.read_excel(pig_folder +'fittedPS/fittedPS_UGM.xlsx')
df_p6 = pd.read_excel(pig_folder +'fittedPS/fittedPS_UGM18.xlsx')


df_list = [df_p1, df_p2, df_p3, df_p4, df_p5, df_p6]
terms_to_remove = ['P1.csv', 'P14.csv', 'P15.csv', 'P18.csv', 'P19.csv', 'P22.csv', 'P23.csv', 'P28.csv', 'P29.csv'] #extraneal dwells we are focusing on physioneal
for i, df in enumerate(df_list):
    print(i)
    df = df[~df.iloc[:,0].isin(terms_to_remove)]
    df_list[i] = df
# df_W0 = df_W0[~df_W0.iloc[:,0].isin(terms_to_remove)]
# df_W1 = df_W1[~df_W1.iloc[:,0].isin(terms_to_remove)]
model_names = ['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM-18']
###############################################
# PLOT OBJECTIVE FN AND COMP TIME
###############################################

# Initialize lists to store mean and std values
mean_obj_fn = []
std_obj_fn = []
mean_comp_time = []
std_comp_time = []

# Compute mean and std for each dataframe
for i, df in enumerate(df_list):
    if i == 0:  # Use C_err for the first model
        mean_obj_fn.append(df['C_err'].mean())
        std_obj_fn.append(df['C_err'].std())
    else:
        mean_obj_fn.append(df['obj_fn'].mean())
        std_obj_fn.append(df['obj_fn'].std())
    mean_comp_time.append(df['comp time'].mean())
    std_comp_time.append(df['comp time'].std())

print(mean_obj_fn, std_obj_fn)
ax1.errorbar(x_values+0.08, mean_obj_fn, yerr=std_obj_fn,  marker='s', color='#CA1F7B', capsize=5, ls = '')
# mean_W0 = df_W0['obj_fn'].mean()
# std_W0 = df_W0['obj_fn'].std()
# ax1.errorbar(1.8, mean_W0, yerr= std_W0, marker='>', color='#8e4585', capsize=5, ls = '' )
# mean_W1 = df_W1['obj_fn'].mean()
# std_W1 = df_W1['obj_fn'].std()
# ax1.errorbar(2.2, mean_W1, yerr= std_W1, marker='<', color='#FF91AF', capsize=5, ls = '' )
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Create a secondary y-axis for 'comp_time'
# ax2 = ax1.twinx()
ax2.errorbar(x_values,np.array(mean_comp_time) / 1000, yerr=np.array(std_comp_time) / 1000, marker='s', color='#CA1F7B', capsize=5, ls = '')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)



#legend
# Create custom legend handles
pig_marker = mlines.Line2D([], [], color='none', marker='s', linestyle='None',
                           markersize=10, label='', markerfacecolor='#CA1F7B', mec = 'None')
human_marker = mlines.Line2D([], [], color='none', marker='o', linestyle='None',
                             markersize=10, label='', markerfacecolor='#0070BB', mec = 'None')

ax2.legend(handles=[pig_marker, human_marker], labels = ['pigs', 'humans'], loc = 'center left', frameon = False)

# # Open the pig image from my computer
# image = open_image_local('pig-head.png')


# # Define the position and size parameters
# image_xaxis = 0.11
# image_yaxis = 0.89
# image_width = 0.05
# image_height = 0.05  # Same as width since our logo is a square

# # Define the position for the image axes
# ax_image = fig.add_axes([image_xaxis,
#                          image_yaxis,
#                          image_width,
#                          image_height]
#                        )

# # Display the image
# ax_image.imshow(image)
# ax_image.axis('off') 

# # Open the pig image from my computer
# image = open_image_local('human-head.png')

# # Define the position and size parameters
# image_xaxis = 0.11
# image_yaxis = 0.84
# image_width = 0.05
# image_height = 0.05  # Same as width since our logo is a square

# # Define the position for the image axes
# ax_image = fig.add_axes([image_xaxis,
#                          image_yaxis,
#                          image_width,
#                          image_height]
#                        )

# # Display the image
# ax_image.imshow(image)
# ax_image.axis('off') 


# Display the plot
plt.tight_layout()
plt.show()
plt.savefig('Figure4.svg', dpi = 600)
