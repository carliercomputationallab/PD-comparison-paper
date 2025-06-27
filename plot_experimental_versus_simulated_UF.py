import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
from fnmatch import fnmatch
import os
import matplotlib.pyplot as plt
import seaborn as sns

"Get all MTACs for the pig in question"
root = 'pigs/patient_files/'
pattern = "*.csv"
patientlist = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            patientlist.append(os.path.join(path, name))
            
df_pigs = pd.read_excel('pigs/fittedPS/fittedPS_TPM.xlsx', index_col = 0)


terms_to_remove = ['P1.csv','P14.csv', 'P15.csv', 'P18.csv', 'P19.csv', 'P22.csv', 'P23.csv', 'P28.csv', 'P29.csv']
# terms_to_remove = ['P28.csv', 'P29.csv']
df_pigs = df_pigs[['Total UF', 'L', 'V_err']]

df_pigs.drop(labels=terms_to_remove, inplace = True)

UF_e = []

mark = []

for pfile in patientlist:
    ind = pfile[33:]  # Adjust index based on file path length
    
    if ind not in terms_to_remove:
        print(ind)
        V = pd.read_csv(pfile,skiprows = range(0,45), delimiter = ",", \
                            encoding= 'unicode_escape')[["IP volume T=0 (mL)","IP volume T=240 (mL)"]].iloc[1]
        if V.isnull().values.any():
            V = pd.read_csv(pfile,skiprows = range(0,45), delimiter = ",", \
                                encoding= 'unicode_escape')[["IP volume T=0 (mL)","IP volume T=240 (mL)"]].iloc[2]
        if V.isnull().values.any():
            V = pd.read_csv(pfile,skiprows = range(0,45), delimiter = ",", \
                                encoding= 'unicode_escape')[["IP volume T=0 (mL)","IP volume T=240 (mL)"]].iloc[0]

        # try:
        #     # Read the csv file, skip rows and handle delimiters
        #     df = pd.read_csv(pfile, skiprows=range(0, 45), delimiter=",", header=0, encoding='unicode_escape')
        #     if 'UFV (mL/h)' in df.columns:
        #         ELAR = float(df['ELAR (mL/h)'].iloc[0])
        #         mark.append('nb')
        #     else:
        #         ELAR = np.nan
        #     # Check for valid UFV values and handle #REF! errors
        #     if 'UFV (mL/h)' in df.columns:
        #         try:
        #             # Attempt to convert to float, catching non-numeric values like '#REF!'
        #             UFV = float(df['UFV (mL/h)'].iloc[2]) * 4  # albumin
        #         except ValueError:  # Catching #REF! or any other invalid data
        #             try:
        #                 UFV = float(df['UFV (mL/h)'].iloc[0]) * 4  # creatinine
        #             except ValueError:
        #                 UFV = np.nan  # If all else fails, set UFV as NaN
        #     else:
        #         try:
        #             UFV = float(df.iloc[25, 7])
        #             ELAR = float(df.iloc[25, 6])/4
        #             mark.append('b')
        #         except ValueError:
        #             UFV = np.nan  # Set UFV as NaN if there's an error
        #             ELAR = np.nan
        #             mark.append('nq')
        # except (ValueError, IndexError) as e:
        #     print(f"Error processing file {pfile}: {e}")
        #     UFV = np.nan  # Set UFV as NaN if there's an error
        #     ELAR = np.nan
        # # print(UFV, ELAR)
        # UF_e.append(UFV-ELAR*4)
        UF_e.append(float(V.iloc[1]) - float(V.iloc[0]))

# Removing NaN values in the collected data
UF_e_clean = [value for value in UF_e if not np.isnan(value)]
df = pd.DataFrame(columns = ['expt', 'sim'])
df['expt'] = UF_e
df['sim'] = df_pigs['Total UF'].values
# df['origin'] = mark
df.index = df_pigs.index
df = df.dropna()
#get human data
datafile = pd.read_excel('human/PD_in_silico_csv_export_/PD_in_silico_export_20240723.xlsx')
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)
patientlist = ['PDPC_059', 'PDPC_236','PDPC_268','PDPC_541','PDPC_572','PDPC_800','PDPC_835','PDPC_890','PDPC_946']




datafile['UF'] =  datafile['Vdrain_t240'] +  datafile['RV240'] -datafile['Vin_t0'] - datafile['RV0'] 

datafile = datafile.loc[patientlist]

df_human  = pd.read_excel('human/fittedPS/fittedPS_TPM.xlsx', index_col = 0)

# Optional: Reorder df_human columns to match df_pigs for consistency
df_human = df_human[['total_ultrafiltration (ml)', 'L', 'V_err']]

df_human.columns = ['Total UF', 'L', 'V_err']
print(df_human['Total UF'])
#%%  
# Set up the figure with a custom layout
fig = plt.figure(figsize=(10, 8))

# Left side subplot (Code 1)
ax_left = plt.subplot2grid((2, 2), (0, 0), rowspan=2)

# Right side subplots (Code 2)
ax_right_top = plt.subplot2grid((2, 2), (0, 1))
ax_right_bottom = plt.subplot2grid((2, 2), (1, 1))

# sns.violinplot(data=[UF_e_clean, df_pigs['Total UF'], datafile['UF'], df_human['Total UF']], inner="point", palette=['#CA1F7B', '#FF91AF', '#0070BB', '#89CFF0'], ax = ax_left)
# ax_left.set_xticks([0, 1, 2, 3], ['Exp', 'Sim','Exp', 'Sim'])
# ax_left.set_ylabel('UF (ml)', fontsize = 14)
# ax_left.tick_params(axis='both', labelsize=14)
# ax_left.spines['top'].set_visible(False)
# ax_left.spines['right'].set_visible(False)
# ax_left.annotate('Pigs', xy=(0.25, -0.1), xycoords='axes fraction', ha='center', fontsize=14, weight='bold')
# ax_left.annotate('Humans', xy=(0.75, -0.1), xycoords='axes fraction', ha='center', fontsize=14, weight='bold')
# ax_left.text(-0.1, 1.05, 'a', transform=ax_left.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')






#%%

# Scatter plot for pigs
ax_left.scatter(df['expt'], df['sim'], c='#CA1F7B', label='pigs')

# Fit a trendline for pigs
pigs_fit = np.polyfit(df['expt'], df['sim'], 1)  
pigs_trendline = np.polyval(pigs_fit, df['expt'])
ax_left.plot(df['expt'], pigs_trendline, c='#CA1F7B', linestyle='-')

# Calculate and display the correlation coefficient for pigs
pigs_corr = np.corrcoef(df['expt'], df['sim'])[0, 1]
ax_left.text(0.05, 0.9, f'Corr (pigs): {pigs_corr:.2f}', transform=ax_left.transAxes, fontsize=12, color='#CA1F7B')

# Scatter plot for humans
ax_left.scatter(datafile['UF'], df_human['Total UF'], c='#0070BB', label='humans')

# Fit a trendline for humans
humans_fit = np.polyfit(datafile['UF'], df_human['Total UF'], 1)  
humans_trendline = np.polyval(humans_fit, datafile['UF'])
ax_left.plot(datafile['UF'], humans_trendline, c='#0070BB', linestyle='-')

# Calculate and display the correlation coefficient for humans
humans_corr = np.corrcoef(datafile['UF'], df_human['Total UF'])[0, 1]
ax_left.text(0.05, 0.85, f'Corr (humans): {humans_corr:.2f}', transform=ax_left.transAxes, fontsize=12, color='#0070BB')

# Labels and styling
ax_left.set_ylabel('Predicted net UF [ml]', fontsize=14)
ax_left.set_xlabel('Measured net UF [ml]', fontsize=14)
ax_left.tick_params(axis='both', labelsize=14)
ax_left.spines['top'].set_visible(False)
ax_left.spines['right'].set_visible(False)
ax_left.text(-0.1, 1.05, 'a', transform=ax_left.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
# ax_left.set_xlim(-500, 1000)
# ax_left.set_ylim(-500, 1000)


#%%
# plt.figure()
# ax_left.scatter(np.arange(len(UF_e_clean)),np.array(UF_e_clean)/df_pigs['Total UF'].to_numpy(), c='#CA1F7B', label = 'pigs')
# ax_left.scatter(np.arange(len(datafile)), datafile['UF'].to_numpy()/ df_human['Total UF'].to_numpy(), c='#0070BB', label = 'humans')
# ax_left.set_ylabel('(Experimental/Predicted) net UF', fontsize = 14)
# ax_left.set_xlabel('pig no.', fontsize = 14)
# ax_left.tick_params(axis='both', labelsize=14)
# ax_left.spines['top'].set_visible(False)
# ax_left.spines['right'].set_visible(False)
# ax_left.legend(frameon = False)
# ax_left.text(-0.1, 1.05, 'a', transform=ax_left.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')


# Create a combined DataFrame for violinplot
df_pigs['Species'] = 'pigs'
df_human['Species'] = 'humans'

# Combine pig and human data
df_combined = pd.concat([df_pigs, df_human], ignore_index=True)

# Plotting for variable L (top right)
sns.violinplot(x='Species', y='L', data=df_combined, inner='point', ax=ax_right_top, palette=['#CA1F7B', '#0070BB'])
ax_right_top.set_ylabel('L [ml/min]', fontsize=14)
ax_right_top.set_xlabel('')
ax_right_top.tick_params(axis='y', labelsize=14)
ax_right_top.spines['top'].set_visible(False)
ax_right_top.spines['right'].set_visible(False)
ax_right_top.text(-0.1, 1.05, 'b', transform=ax_right_top.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

# Plotting for variable V_err (bottom right)
sns.violinplot(x='Species', y='V_err', data=df_combined, inner='point', ax=ax_right_bottom, palette=['#CA1F7B', '#0070BB'])
ax_right_bottom.set_ylabel('Volume error [-]', fontsize=14)
ax_right_bottom.set_xlabel('')
ax_right_bottom.tick_params(axis='y', labelsize=14)
ax_right_bottom.spines['top'].set_visible(False)
ax_right_bottom.spines['right'].set_visible(False)
ax_right_bottom.text(-0.1, 1.09, 'c', transform=ax_right_bottom.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.show()

plt.savefig('Figure7.svg', dpi = 600)