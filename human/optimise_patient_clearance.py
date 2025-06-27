# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:14:18 2024
try to find which treatment option gives the best "treatment score"
@author: P70073624
"""



# Importing necessary libraries and modules

import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns

st = time.time()
#%%
# Defining a function named "objective"
# This function is used as an objective function to be minimized in optimization
def objective(x, predicted_cd, Cp, V,  Vr, V_fill,  LS, t):
    '''The objective function needed to be minimised'''
    # ...
    # The function has various calculations and returns multiple values
    # It calculates and returns the sum of squared differences between two sets of data
    # It also returns predicted values and other variables

    V = np.zeros(t+1)
    V[0] = V_fill+Vr
    predicted_cd, UF, V = rk(t, x, predicted_cd, Cp, V,  Vr, V_fill, LS)
    
    
    return ( predicted_cd, UF, V)
    
    

#%% OBJECTIVE FUNCTION FOR MINIMISE CANNOT HAVE MORE THAN ONE RETURN OUTPUT, so we 

#%%
#Runge-Kutta
def rk(t, x, predicted_cd, Cp, V,  Vr, V_fill, LS):
    # print(x)
    Jv = []
    
    
    for timestep in range(0,t): 
        # print(predicted_V[timestep])
        cd = predicted_cd.loc[timestep]
        
        "Apply Runge Kutta Formulas to find next value of y"
        k1, v1 = comdxdt(cd, timestep, x,  predicted_cd, Cp, V, Vr, V_fill, LS)
        k2, v2 = comdxdt(cd + 0.5  *k1, timestep, x,  predicted_cd, Cp, V, Vr, V_fill, LS)
        k3, v3 = comdxdt(cd + 0.5  *k2, timestep, x,  predicted_cd, Cp, V,  Vr, V_fill, LS)
        k4, v4 = comdxdt(cd + k3, timestep, x,  predicted_cd, Cp, V,  Vr, V_fill, LS)
        
        # Update next value of y
        cd = cd + (1 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        Jv.append( (1 / 6.0)*(v1 + 2 * v2 + 2 * v3 + v4))
        #print(UF)
        predicted_cd.loc[timestep+1] = cd
        V[timestep+1] = V[timestep] + (1 / 6.0)*(v1 + 2 * v2 + 2 * v3 + v4)
    # print(x, predicted_V[t])            
    return predicted_cd, Jv, V

# the differential equations
def comdxdt(cd, t, x, predicted_cd, Cp, V,  Vr, V_fill, LS):
    
    solutes = ["Urea", "Creatinine", "Glucose"]
    
    #MTAC for small pores
    PS_s = x[0:3] * 0.998 * abya0_s #fraction small pore surface area - Rippe, A THREE-PORE MODEL OF PERITONEAL TRANSPORT table 1

    #MTAC for large pores
    PS_l = x[0:3] * 0.002 * abya0_l# Ref: two pore model-oberg,rippe, table 1, A0L/A0 value
    
    # alpha[0] = x[3]
    # alpha[1] = 1- x[3] - alpha[2]
    
    
    #fraction of peritoneal membrane in contact with the dialysis fluid
    af = 16.18 * (1 - np.exp (-0.77*V[t]))/13.3187
    
    #hydrostatic pressure difference
    delP = delP0 - ((V[t] - (V_fill+Vr))/490)
    
    #peritoneal concentration gradient
    pr = [phi[i]*RT * (Cp[i]-cd[i]) for i in range(len(solutes))]
    sigmas_pr = sum([phi[i] * sigma_s[i]*pr[i] for i in range(len(solutes))])
    sigmal_pr = sum([phi[i] * sigma_l[i]*pr[i] for i in range(len(solutes))])

    # #volumetric flows across the pores
    J_vC = af*alpha[0]*LS*(delP - sum(pr)) #ml/min
    J_vS = af*alpha[1]*LS*(delP  - sigmas_pr) #ml/min
    J_vL = af*alpha[2]*LS*(delP - sigmal_pr) #ml/min

    # #Peclet numbers
    Pe_s = np.array([J_vS  * (1 - sigma_s[i])/(af*PS_s[i]) for i in range(len(solutes))])
    Pe_s[np.isinf(Pe_s)] = 0
    Pe_l = np.array([J_vL  * (1 - sigma_l[i])/(af*PS_l[i]) for i in range(len(solutes))])
    Pe_l[np.isinf(Pe_l)] = 0
    
    
    # #solute flow rate
    J_sS = (J_vS*(1-sigma_s)*(Cp-cd*np.exp(-Pe_s))/(1-np.exp(-Pe_s))).ravel()
    J_sL = (J_vL*(1-sigma_l)*(Cp-cd*np.exp(-Pe_l))/(1-np.exp(-Pe_l))).ravel()
    
    dvdt = J_vC+J_vS+J_vL-L

    dxdt = ((J_sS + J_sL)/V[t]-np.array(cd)*(J_vC + J_vS + J_vL-L)/V[t]).ravel()

    return (dxdt, dvdt)


    
#%%
# get patient specific data
datafile = pd.read_csv('PD_in_silico_csv_export_/PD_in_silico_export_20240723.csv', sep = ";")
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)
patientlist = datafile.index

data = {}

for i,ind in enumerate(datafile.index):
    
    
    Vfill = datafile.loc[ind, 'Vin_t0']
    
    RV = 200
    
    V_drain = datafile.loc[ind, 'Vdrain_t240']
        
    
    
    cp = datafile.loc[ind, ['Ur_t2_B','Cr_t2_B', 'Gluc_t2_B']]
    cp.index = ["Urea", "Creatinine", "Glucose"]
    if cp.loc['Glucose'] < 0:
        cp.loc['Glucose'] = 6.0
        
    cp.loc['Creatinine'] *=0.001
    
    df_V = np.array(datafile.loc[ind,['Vin_t0', 'Vdrain_t240']]) #mL
    df_V[1] += RV #assuming a residual volume of 200 ml
    #Linear interpolation to find values of V at all times
    f_V = interp1d([0,240], df_V)
    interpolated_V = f_V(range(0,241))

    V = interpolated_V

    Vr = RV  #mL
    
    # LS = LPS.loc[ind, 'LpS'] #ml/min mmHg
    
    # if np.isnan(LS):
    #     LS = 0.02
        
    # LPS.loc[ind, 'LpS'] = LS
    
    LS = 0.011
    
    pat_d = { 'cp': cp, 'V': V, 'Vfill': Vfill, 'V_drain': V_drain, 'Vr': Vr, 'LS':LS}
    
    data[ind] = pat_d
  
# ultrafiltration coefficient


# fractional pore coefficients
alpha = [0.020, 0.900, 0.080]

# Initital hydrostatic pressure difference
delP0 = 8 #mmHg

# constant
RT = 19.3 #mmHg per mmol/l
  
#small pore radius
rs = 43 # Angstrom
#large pore radius
rl = 250

solutes = ["Urea", "Creatinine", "Glucose"]
#radius of molecules
r = np.array([ 2.6, 3.0, 3.7]) #the phosphate radius is approximated from its topological surface area
#for radius, new paper by Oberg - https://journals.sagepub.com/doi/suppl/10.1177/08968608211069232
#constants to calculate sigma
gamma_s = r/rs
gamma_l = r/rl

# lymphatic absorption rate
L = 0.3 #ml/min 

abya0_s = 1+9/8*gamma_s*np.log(gamma_s)-1.56034*gamma_s+0.528155*gamma_s**2+\
    1.91521*gamma_s**3-2.81903*gamma_s**4+0.270788*gamma_s**5+1.10115*gamma_s**6+ 0.435933*gamma_s**7 #eq 21 two pore Ficoll
abya0_l = 1+9/8*gamma_l*np.log(gamma_l)-1.56034*gamma_l+0.528155*gamma_l**2+\
    1.91521*gamma_l**3-2.81903*gamma_l**4+0.270788*gamma_l**5+1.10115*gamma_l**6+ 0.435933*gamma_l**7

#Osmotic reflection coefficients
sigma_s = np.zeros(len(solutes))
sigma_l = np.zeros(len(solutes))
sigma = np.zeros(len(solutes))

for i in range(len(solutes)):
    sigma_s[i] = 16/3 * (gamma_s[i])**2 - 20/3 * (gamma_s[i])**3 + 7/3 * (gamma_s[i])**4
    sigma_l[i] = 16/3 * (gamma_l[i])**2 - 20/3 * (gamma_l[i])**3 + 7/3 * (gamma_l[i])**4
    sigma[i] = alpha[0] + alpha[1] * sigma_s[i] + alpha[2] * sigma_l[i]

# dissociation factor
phi = np.array([1, 1,  1])





MTAC = pd.read_excel('fittedPS/fittedPS_TPM.xlsx', index_col = 0)    
for ind in ['754']:

    pfile = f'PDPC_{ind}'
    Cp = data[pfile]['cp']
    V = data[pfile]['V']
    Vr = 200
    V_fill = data[pfile]['Vfill']
    LS = data[pfile]['LS']
    V_drain = data[pfile]['V_drain']
    
    predicted_cd = pd.DataFrame(columns= ["Urea", "Creatinine", "Glucose"], dtype = float)
    
    predicted_cd.loc[0]= np.array([0, 0, 0])
    
    result = pd.DataFrame(columns = ['t_final', 'glu_0', 'urea', 'crea', 'UF', 'glucose_abs'])

    # for 3 dwell times and 3 glucose concentraiton run TPM to get UF and clearance of urea and creatinine
    for j, t_final in enumerate([60, 120, 480]):
        for i, glu in enumerate([75, 126, 214]):
            predicted_cd.loc[0]= np.array([0, 0, glu])
            x = np.array(MTAC.loc[pfile].iloc[0:3])
            predicted_cd, UF, V = objective(x, predicted_cd, Cp, V,  Vr, V_fill,  LS, t_final)
            urea = (predicted_cd.loc[t_final, 'Urea']* V[-1]- predicted_cd.loc[0,'Urea']* (V_fill+Vr))/Cp['Urea']
            crea = (predicted_cd.loc[t_final, 'Creatinine']* V[-1]- predicted_cd.loc[0,'Creatinine']*(V_fill+Vr))/Cp['Creatinine']
            glucose_abs = predicted_cd.loc[0,'Glucose']* (V_fill+Vr)-predicted_cd.loc[t_final, 'Glucose']* V[-1]
            urea /= glucose_abs
            
            uf = sum(UF)/t_final/glucose_abs
            result.loc[3*j+i] = [t_final, glu, urea, crea, uf, glucose_abs]
  
    #normalise the data so obtained between 0 and 1 (scaled per column of urea crea and UF)
    scaler = MinMaxScaler() 
    arr_scaled = scaler.fit_transform(result) 
    df_scaled = pd.DataFrame(arr_scaled, columns=result.columns,index=result.index)
    result['urea_scaled'] = df_scaled['urea']
    result['crea_scaled'] = df_scaled['crea']
    result['UF_scaled'] = df_scaled['UF']
    result['glu_abs_scaled'] = df_scaled['glucose_abs']
    
    # rearrange dfs
    pivot_df1 = pd.DataFrame(columns=[75, 126, 214], index=[60, 120, 240])
    pivot_df2 = pd.DataFrame(columns=[75, 126, 214], index=[60, 120, 240])
    pivot_df3 = pd.DataFrame(columns=[75, 126, 214], index=[60, 120, 240])
    
    pivot_df = [pivot_df1, pivot_df2, pivot_df3]
    
    # patients can have different priority, this i am parameterising as weights
    # multiplying weights withe the normalised values we can calculate treatment score
    # treatmetn score = weight clearance * urea clearance + weight UF * UF
    weights = [0.3, 0.5, 0.7]
    for i, wt in enumerate(weights):
        wt_clearance = wt
        wt_UF = 1 - wt_clearance
        result['Treatment_function'] = wt_clearance * result['urea_scaled'] + wt_UF * result['UF_scaled']
    
        # Pivot the DataFrame to get it into the shape we need for the heatmap
        # We want `t_final` as rows, `glu_0` as columns, and `Treatment_function` as the values
        pivot_df[i] = result.pivot("t_final", "glu_0", "Treatment_function")    

    # Finding global min and max values across all datasets
    global_min = min(pivot_df[0].min().min(), pivot_df[1].min().min(), pivot_df[2].min().min())
    global_max = max(pivot_df[0].max().max(), pivot_df[1].max().max(), pivot_df[2].max().max())

    
    # Plotting the heatmaps with a unified color scale
    plt.figure(2)
    fig, axs = plt.subplots(1, 3, figsize=(14, 6), sharey=True, sharex=True)
    cbar_ax = fig.add_axes([.91, .3, .03, .4])  # Positioning the color bar
    tick_positions = [0.5, 1.5, 2.5]  # Position of ticks corresponding to the columns in the dataframe
    tick_labels = [1.36, 2.27, 3.86]  # Custom labels for the ticks
    annot_font_size = 16
    index = ['(a)', '(b)', '(c)']
    for i, ax in enumerate(axs):
        sns.heatmap(pivot_df[i], annot=True, cmap="viridis", fmt=".2f", ax=ax,
                    cbar=i == 0, cbar_ax=None if i else cbar_ax, vmin=global_min, vmax=global_max,
                    annot_kws={"size": annot_font_size})
        ax.set_title(f'{index[i]} $w_{{\mathrm{{clearance}}}} = {weights[i]}$, $w_{{\mathrm{{UF}}}} = {1 - weights[i]:.1f}$', fontsize=16)
        ax.set_xticks(tick_positions)  # Setting the positions for the ticks
        ax.set_xticklabels(tick_labels, fontsize=16)
        ax.set_xlabel('')
        ax.set_ylabel('')
    axs[0].set_yticks(tick_positions)
    axs[0].set_yticklabels([1, 2, 4], fontsize=16)
    fig.supxlabel('Initial glucose (% w/v)', fontsize=16)
    fig.supylabel('Dwell time (hrs)', fontsize=16)
    
    plt.tight_layout(rect=[0, 0, .9, 1])  # Adjust the layout to make room for the colorbar
    plt.subplots_adjust(top=0.912,
                        bottom=0.172,
                        left=0.084,
                        right=0.889,
                        hspace=0.2,
                        wspace=0.154)
    plt.show()
    plt.savefig(f'final-plots/personalised_treatment_PDPC_{ind}.png', dpi = 600)

    