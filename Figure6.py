import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Sample data for pig dataframe
models = ['TPM', 'GM', 'WM', 'SWM', 'UGM', 'UGM-18']
df_pig = pd.read_excel('pigs/population_average_pigs.xlsx', index_col = 0)
df_pig = df_pig[models]

df_human = pd.read_excel('human/population_average_human.xlsx', index_col = 0)
df_human = df_human[models]
# Create a new dataframe by alternating rows from pig and human dataframes
new_data = {col: [] for col in df_pig.columns}
new_index = []

for idx in df_pig.index:
    # Add row from pig dataframe
    for col in df_pig.columns:
        new_data[col].append(df_pig.at[idx, col])
    new_index.append(idx)
    
    # Add corresponding row from human dataframe, if exists
    if idx in df_human.index:
        for col in df_pig.columns:
            new_data[col].append(df_human.at[idx, col])
        new_index.append(idx + '_human')
    else:
        for col in df_pig.columns:
            new_data[col].append(np.nan)
        new_index.append(idx + '_human')

# Create the final dataframe
df_combined = pd.DataFrame(new_data, index=new_index)

# save the combined dataframe
df_combined.to_excel('population_average_pigs_human.xlsx')



# Function to extract mean and std dev
def extract_mean_std(data_str):
    if isinstance(data_str, (float, int)):
        if np.isnan(data_str):
            return 0, 0
        else:
            return data_str, 0
    elif data_str:
        mean, std = data_str.split(' ± ')
        return float(mean), float(std)
    elif np.isnan(data_str):
        return 0, 0
    else:
        return 0, 0
solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
# Plotting MTAC
# Create a 3x2 plot
Literature = [15.11, 5.682, 1.936, 4.171, 7.319, 24.352]
fig, axes = plt.subplots(3, 2, figsize=(12, 12), sharex = True)
axes = axes.flatten()
index = df_combined.index[:12]
for i, (ax, idx) in enumerate(zip(axes, range(0, len(index), 2))):
    row1, row2 = index[idx], index[idx+1]
    means_pigs, stds_pigs, means_humans, stds_humans = [], [], [], []
    x_pigs = np.arange(-0.25, 0.3, 0.1) # All Pigs bars around x=0
    x_humans = np.arange(0.75, 1.3, 0.1)  # All Humans bars around x=1
    width = 0.1  # Width of the bars
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    for j, model in enumerate(df_combined.keys()):
        mean1, std1 = extract_mean_std(df_combined[model][idx])
        mean2, std2 = extract_mean_std(df_combined[model][idx+1])
        # Plot Pigs df_combined
        ax.bar(x_pigs[j], mean1, width, yerr=std1, capsize=5, align='center', color = colors[j], edgecolor = 'k', label = model)
        # Plot Humans df_combined
        ax.bar(x_humans[j], mean2, width, yerr=std2, capsize=5, align='center', color = colors[j], edgecolor = 'k')
        # Add literature value lines and ±10% lines
        lit_value = Literature[i]
        ax.hlines(lit_value, 0.5, 1.5, colors='gray', linestyles='dashed', linewidth=1.5)
        ax.hlines(lit_value * 1.1, 0.5, 1.5, colors='gray', linestyles='dotted', linewidth=1)
        ax.hlines(lit_value * 0.9, 0.5, 1.5, colors='gray', linestyles='dotted', linewidth=1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Set x-ticks and labels
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Pigs', 'Humans'], fontsize = 14)
    ax.set_title(solutes[i])
    ax.tick_params('y', labelsize = 14)
axes[2].set_ylim(-12.5, 6.8)
axes[3].legend(frameon = False)
fig.supylabel('MTAC [ml/min]', fontsize = 14)
plt.tight_layout()
plt.show()
plt.savefig('Figure6.svg', dpi = 600)

# Plotting fct
# Create a 3x2 plot
fig, axes = plt.subplots(3, 2, figsize=(12, 12), sharex = True)
axes = axes.flatten()
index = df_combined.index[12:24]
for i, (ax, idx) in enumerate(zip(axes, range(0, len(index), 2))):
    row1, row2 = index[idx], index[idx+1]
    means_pigs, stds_pigs, means_humans, stds_humans = [], [], [], []
    x_pigs = np.arange(-0.25, 0.3, 0.1) # All Pigs bars around x=0
    x_humans = np.arange(0.75, 1.3, 0.1)  # All Humans bars around x=1
    width = 0.1  # Width of the bars
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    for j, model in enumerate(df_combined.keys()):
        mean1, std1 = extract_mean_std(df_combined[model][row1])
        mean2, std2 = extract_mean_std(df_combined[model][row2])
        # Plot Pigs df_combined
        ax.bar(x_pigs[j], mean1, width, yerr=std1, capsize=5, align='center', color = colors[j], edgecolor = 'k', label = model)
        # Plot Humans df_combined
        ax.bar(x_humans[j], mean2, width, yerr=std2, capsize=5, align='center', color = colors[j], edgecolor = 'k')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Set x-ticks and labels
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Pigs', 'Humans'], fontsize = 14)
    ax.set_title(solutes[i])
    ax.tick_params('y', labelsize = 14)
    
axes[3].legend(frameon = False)
fig.supylabel('fct', fontsize = 14)
plt.tight_layout()
plt.show()
plt.savefig('population_average_fct.png', dpi = 600)

# Plotting sico
# Create a 3x2 plot
fig, axes = plt.subplots(3, 2, figsize=(12, 12), sharex = True)
axes = axes.flatten()
index = df_combined.index[24:36]
for i, (ax, idx) in enumerate(zip(axes, range(0, len(index), 2))):
    row1, row2 = index[idx], index[idx+1]
    means_pigs, stds_pigs, means_humans, stds_humans = [], [], [], []
    x_pigs = np.arange(-0.25, 0.3, 0.1) # All Pigs bars around x=0
    x_humans = np.arange(0.75, 1.3, 0.1)  # All Humans bars around x=1
    width = 0.1  # Width of the bars
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    for j, model in enumerate(df_combined.keys()):
        mean1, std1 = extract_mean_std(df_combined[model][row1])
        mean2, std2 = extract_mean_std(df_combined[model][row2])
        # Plot Pigs df_combined
        ax.bar(x_pigs[j], mean1, width, yerr=std1, capsize=5, align='center', color = colors[j], edgecolor = 'k', label = model)
        # Plot Humans df_combined
        ax.bar(x_humans[j], mean2, width, yerr=std2, capsize=5, align='center', color = colors[j], edgecolor = 'k')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Set x-ticks and labels
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Pigs', 'Humans'], fontsize = 14)
    ax.set_title(solutes[i])
    ax.tick_params('y', labelsize = 14)
axes[3].legend(frameon = False)
fig.supylabel('SiCo', fontsize = 14)
plt.tight_layout()
plt.show()
plt.savefig('population_average_sico.png', dpi = 600)