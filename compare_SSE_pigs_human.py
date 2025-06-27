# -s- coding: utf-8 -s-
"""
Created on Tue Apr  9 10:07:48 2024

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

####################################################
#############              PIGS          ###########
####################################################

# Read the specific solute error collected from sse python files
df_p1 = pd.read_excel('pigs/fittedPS/sse_TPM.xlsx', index_col = 0)
df_p2 = pd.read_excel('pigs/fittedPS/sse_UGM.xlsx', index_col = 0)
df_p3 = pd.read_excel('pigs/fittedPS/sse_UGM18.xlsx', index_col = 0)
df_p4 = pd.read_excel('pigs/fittedPS/sse_GM.xlsx', index_col = 0)
df_p5 = pd.read_excel('pigs/fittedPS/sse_WM_F0p5.xlsx', index_col = 0)
df_W0 = pd.read_excel('pigs/fittedPS/sse_WM_F0.xlsx', index_col = 0)
df_W1 = pd.read_excel('pigs/fittedPS/sse_WM_F1.xlsx', index_col = 0)
df_p6 = pd.read_excel('pigs/fittedPS/sse_SWM.xlsx', index_col = 0)
labels = ['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM-18']

dataframes = [df_p1, df_p4, df_p5, df_p6, df_p2, df_p3]
terms_to_remove = ['P1.csv', 'P14.csv', 'P15.csv', 'P18.csv', 'P19.csv', 'P22.csv', 'P23.csv', 'P28.csv', 'P29.csv'] #extraneal dwells we are focusing on physioneal
for i, df in enumerate(dataframes):
    dataframes[i] = df[~df.index.isin(terms_to_remove)]
df_W0 = df_W0[~df_W0.index.isin(terms_to_remove)]
df_W1 = df_W1[~df_W1.index.isin(terms_to_remove)]
# Plotting pigs
solute_columns = df_p1.columns
num_rows = 2
num_cols = 3

fig, axs = plt.subplots(num_rows, num_cols, figsize=(9, 6), sharex = True)
print('PIGS')
for i, solute in enumerate(solute_columns):
    row_index = i // num_cols
    col_index = i % num_cols
    
    means = [df[solute].mean() for df in dataframes]
    stds = [df[solute].std() for df in dataframes]
    
    print(solute, f'{means} ± {stds}')
    x_pos = np.arange(len(labels))

    axs[row_index, col_index].errorbar(x_pos-0.1, means, yerr=stds, color = '#CA1F7B', marker = 's', ecolor='#CA1F7B', capsize=5, ls = '')
    axs[row_index, col_index].set_title(solute, fontsize =12)
    # axs[row_index, col_index].errorbar(1.8, df_W0[solute].mean(), df_W0[solute].std(),marker='>', color='#CA1F7B', capsize=5, ls = '')
    # axs[row_index, col_index].errorbar(2.2, df_W1[solute].mean(), df_W1[solute].std(),marker='<', color='#8e4585', capsize=5, ls = '')
    
    # comment out to plot the individual data points
    # for j, df in enumerate(dataframes):
    #     axs[row_index, col_index].scatter([x_pos[j]] s len(df), df[solute], s=10, alpha=0.5, color='gray')



####################################################
#############             HUMAN          ###########
####################################################

# Read the specific solute error collected from sse python files
df_h1 = pd.read_excel('human/fittedPS/sse_TPM.xlsx', index_col = 0)
df_h2 = pd.read_excel('human/fittedPS/sse_UGM.xlsx', index_col = 0)
df_h3 = pd.read_excel('human/fittedPS/sse_UGM18.xlsx', index_col = 0)
df_h4 = pd.read_excel('human/fittedPS/sse_GM.xlsx', index_col = 0)
df_h5 = pd.read_excel('human/fittedPS/sse_WM_F0p5.xlsx', index_col = 0)
# df_W0 = pd.read_excel('human/fittedPS/sse_WM_F0.xlsx', index_col = 0)
# df_W1 = pd.read_excel('human/fittedPS/sse_WM_F1.xlsx', index_col = 0)
df_h6 = pd.read_excel('human/fittedPS/sse_SWM.xlsx', index_col = 0)
labels = ['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM-18']
#insert others here
dataframes = [df_h1, df_h4, df_h5, df_h6, df_h2, df_h3]

# Plotting
solute_human = df_h1.columns


print('HUMAN')

for i, solute in enumerate(solute_columns):
    row_index = i // num_cols
    col_index = i % num_cols
    if solute in solute_human:
        means = [df[solute].mean() for df in dataframes]
        stds = [df[solute].std() for df in dataframes]
        
        print(solute, f'{means} ± {stds}')
        x_pos = np.arange(len(labels))
    
        axs[row_index, col_index].errorbar(x_pos+0.1, means, yerr=stds, color = '#0070BB', marker = 'o',  ecolor='#0070BB', capsize=5, ls = '')
        # axs[row_index, col_index].errorbar(1.8, df_W0[solute].mean(), df_W0[solute].std(),marker='>', color='#103783', capsize=5, ls = '')
        # axs[row_index, col_index].errorbar(2.2, df_W1[solute].mean(), df_W1[solute].std(),marker='<', color='#432371', capsize=5, ls = '')
        
    
for ax in axs.flat:
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation = 30)
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax.tick_params(axis = 'both', labelsize = 12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# Create custom legend handles
pig_marker = mlines.Line2D([], [], color='none', marker='s', linestyle='None',
                            markersize=12, label='', markerfacecolor='#CA1F7B', mec = 'None')
human_marker = mlines.Line2D([], [], color='none', marker='o', linestyle='None',
                              markersize=12, label='', markerfacecolor='#0070BB', mec = 'None')

axs[0,0].legend(handles=[pig_marker, human_marker], labels = ['pigs', 'humans'], loc = 'upper left', frameon = False)

# Open the pig image from my computer
# image = open_image_local('pig-head.png')


# # Define the position and size parameters
# image_xaxis = 0.14
# image_yaxis = 0.89
# image_width = 0.05
# image_height = 0.05  # Same as width since our logo is a square

# # Define the position for the image axes
# ax_image = fig.add_axes([image_xaxis,
#                           image_yaxis,
#                           image_width,
#                           image_height]
#                         )

# # Display the image
# ax_image.imshow(image)
# ax_image.axis('off') 

# # Open the pig image from my computer
# image = open_image_local('human-head.png')

# # Define the position and size parameters
# image_xaxis = 0.14
# image_yaxis = 0.85
# image_width = 0.045
# image_height = 0.045  # Same as width since our logo is a square

# # Define the position for the image axes
# ax_image = fig.add_axes([image_xaxis,
#                           image_yaxis,
#                           image_width,
#                           image_height]
#                         )

# # Display the image
# ax_image.imshow(image)
# ax_image.axis('off') 
fig.supylabel('SSE', fontsize = 12)
plt.tight_layout()
plt.show()

plt.savefig('Figure5.svg', dpi = 600)


