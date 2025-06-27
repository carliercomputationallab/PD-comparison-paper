# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:24:46 2024

@author: P70073624
"""



# Importing necessary libraries and modules

import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


st = time.time()
#%%
# Defining a function named "objective"
# This function is used as an objective function to be minimized in optimization
def objective(x, predicted_cd, predicted_V, Cp, V, df_cd, LS):
    '''The objective function needed to be minimised'''
    # ...
    # The function has various calculations and returns multiple values
    # It calculates and returns the sum of squared differences between two sets of data
    # It also returns predicted values and other variables
    
    
    t = 240

    predicted_cd, UF = rk(t, x, predicted_cd, predicted_V, Cp, df_cd, LS, V[0])
    df2 = ((df_cd/Cp-predicted_cd.loc[df_cd.index]/Cp)**2).astype('float64')
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
    
    return np.sqrt(df2.sum()/len(df_cd)), predicted_cd
    

#%%
#Runge-Kutta
def rk(t, x, predicted_cd, predicted_V, Cp, df_cd, LS, V0):
    # print(x)
    Jv = []
    for timestep in range(0,t): 
        # print(predicted_V[timestep])
        cd = np.array(predicted_cd.loc[timestep])
        Vp = predicted_V[timestep]
        #we have predicted_V as an array
        #since we start from 20 it doesnt work with the array unless we subtract 20
        "Apply Runge Kutta Formulas to find next value of y"
        k1, v1 = comdxdt(cd, timestep, x,  predicted_cd, Cp, Vp, df_cd, LS, V0)
        k2, v2 = comdxdt(cd + 0.5  *k1, timestep, x,  predicted_cd, Cp, Vp + 0.5*v1, df_cd, LS, V0)
        k3, v3 = comdxdt(cd + 0.5  *k2, timestep, x,  predicted_cd, Cp, Vp + 0.5*v2, df_cd, LS, V0)
        k4, v4 = comdxdt(cd + k3, timestep, x,  predicted_cd, Cp, Vp + v3, df_cd, LS, V0)
        
        # Update next value of y
        cd = cd + (1 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        Jv.append( (1 / 6.0)*(v1 + 2 * v2 + 2 * v3 + v4))
        #print(UF)
        predicted_cd.loc[timestep+1] = cd
        predicted_V[timestep+1] = predicted_V[timestep] + (1 / 6.0)*(v1 + 2 * v2 + 2 * v3 + v4) 
    # print(x, predicted_V[t])            
    return predicted_cd, Jv

# the differential equations
def comdxdt(cd, t, x, predicted_cd, Cp, Vp, df_cd, LS, V0):
    
    solutes = ["Urea", "Creatinine", "Glucose", "Sodium"]
    
    #MTAC for small pores
    PS_s = x[0:4] * 0.998 * abya0_s #fraction small pore surface area - Rippe, A THREE-PORE MODEL OF PERITONEAL TRANSPORT table 1

    #MTAC for large pores
    PS_l = x[0:4] * 0.002 * abya0_l# Ref: two pore model-oberg,rippe, table 1, A0L/A0 value
    
    # alpha[0] = x[4]
    # alpha[1] = 0.92-x[4]
    # alpha[2] = 0.08
    
    L = x[4]
    
    
    #fraction of peritoneal membrane in contact with the dialysis fluid
    af = 16.18 * (1 - np.exp (-0.00077*Vp))/13.3187
    
    #hydrostatic pressure difference
    delP = delP0 - ((Vp - V0)/490)
    
    #peritoneal concentration gradient
    pr = [phi[i]*RT * (Cp[i]-cd[i]) for i in range(len(solutes))]
    sigmas_pr = sum([phi[i] * sigma_s[i]*pr[i] for i in range(len(solutes))])
    sigmal_pr = sum([phi[i] * sigma_l[i]*pr[i] for i in range(len(solutes))])

    # #volumetric flows across the pores
    J_vC = af*alpha[0]*LS*(delP - sum(pr) - 22) #ml/min
    J_vS = af*alpha[1]*LS*(delP  - sigmas_pr-0.963*22) #ml/min
    J_vL = af*alpha[2]*LS*(delP - sigmal_pr-0.083*22) #ml/min

    # #Peclet numbers
    Pe_s = np.array([J_vS  * (1 - sigma_s[i])/(af*PS_s[i]) for i in range(len(solutes))])
    Pe_s[np.isinf(Pe_s)] = 0
    Pe_l = np.array([J_vL  * (1 - sigma_l[i])/(af*PS_l[i]) for i in range(len(solutes))])
    Pe_l[np.isinf(Pe_l)] = 0

    
    # #solute flow rate
    J_sS = (J_vS*(1-sigma_s)*(Cp-cd*np.exp(-Pe_s))/(1-np.exp(-Pe_s))).ravel()
    J_sL = (J_vL*(1-sigma_l)*(Cp-cd*np.exp(-Pe_l))/(1-np.exp(-Pe_l))).ravel()
    
    dvdt = J_vC+J_vS+J_vL-L

    dxdt = ((J_sS + J_sL)/Vp-np.array(cd)*(J_vC + J_vS + J_vL)/Vp).ravel()

    return (dxdt, dvdt)



    
#%%
datafile = pd.read_excel('PD_in_silico_csv_export_/PD_in_silico_export_20240723.xlsx')
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)

patientlist = ['PDPC_059', 'PDPC_236','PDPC_268','PDPC_541','PDPC_572','PDPC_800','PDPC_835','PDPC_890','PDPC_946']

#create a dictionary based on the patient file provided
data = {}

for i,ind in enumerate(patientlist):
    
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
        
    RV0 = datafile.loc[ind, "RV0"]
    RV240 = datafile.loc[ind, "RV240"]
    
    ur_0 = [datafile.loc[ind, 'Ur_night']*RV0/(Vfill+RV0) if datafile.loc[ind, 'Ur_night'] > 0 else 0]
    cr_0 = [datafile.loc[ind, 'Cr_night']*RV0/(Vfill+RV0)*0.001 if datafile.loc[ind, 'Cr_night'] > 0 else 0]
    glu_0 = [214 if datafile.loc[ind, 'Gluc_PET'] == 3.86 else  126 if datafile.loc[ind, 'Gluc_PET'] == 2.27 else  75 ]
    na_0 = [132]
    df_cd.loc[0] = np.concatenate((ur_0, cr_0, glu_0, na_0))
    
    #blood plasma concentration from datafiles
    cp = datafile.loc[ind, ['Ur_t2_B','Cr_t2_B', 'Gluc_t2_B', 'Na_t2_B']]
    cp.index = ["Urea", "Creatinine", "Glucose", "Sodium"]
    if cp.loc['Glucose'] < 0:
        cp.loc['Glucose'] = 6.0
        
    cp.loc['Creatinine'] *=0.001
    cp.loc['Sodium'] *=0.94
    cp = np.array(cp)*1/(0.984 - 0.000718*70) #70 g/L is mean human total protein blood levels
    
    # Volume from datafiles
    df_V = np.array(datafile.loc[ind,['Vin_t0', 'Vdrain_t240']]) #mL
    df_V[0] += RV0 #assuming a residual volume of 200 ml
    df_V[1] += RV240 #assuming a residual volume of 200 ml
    #Linear interpolation to find values of V at all times
    f_V = interp1d([0,240], df_V)
    interpolated_V = f_V(range(0,241))

    V = interpolated_V

    # Set the ultrafiltraiton coefficients
    LS = 0.074 #ml/min mmHg

    pat_d = {'df_cd': df_cd, 'cp': cp, 'V': V,  'LS':LS}
    
    data[ind] = pat_d
  

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

solutes = ["Urea", "Creatinine", "Glucose", "Sodium"]
#radius of molecules
r = np.array([ 2.6, 3.0, 3.7, 2.3]) #the phosphate radius is approximated from its topological surface area
#for radius, new paper by Oberg - https://journals.sagepub.com/doi/suppl/10.1177/08968608211069232
#constants to calculate sigma
gamma_s = r/rs
gamma_l = r/rl

# lymphatic absorption rate
L = 0.3 #ml/min 

abya0_s = (1-gamma_s)**4.5/(1-0.3956*gamma_s+1.0616*gamma_s**2)
abya0_l = (1-gamma_l)**4.5/(1-0.3956*gamma_l+1.0616*gamma_l**2)

# abya0_s = 1+9/8*gamma_s*np.log(gamma_s)-1.56034*gamma_s+0.528155*gamma_s**2+\
#     1.91521*gamma_s**3-2.81903*gamma_s**4+0.270788*gamma_s**5+1.10115*gamma_s**6+ 0.435933*gamma_s**7 #eq 21 two pore Ficoll
# abya0_l = 1+9/8*gamma_l*np.log(gamma_l)-1.56034*gamma_l+0.528155*gamma_l**2+\
#     1.91521*gamma_l**3-2.81903*gamma_l**4+0.270788*gamma_l**5+1.10115*gamma_l**6+ 0.435933*gamma_l**7

#Osmotic reflection coefficients
sigma_s = np.zeros(len(solutes))
sigma_l = np.zeros(len(solutes))
sigma = np.zeros(len(solutes))

for i in range(len(solutes)):
    sigma_s[i] = 16/3 * (gamma_s[i])**2 - 20/3 * (gamma_s[i])**3 + 7/3 * (gamma_s[i])**4
    sigma_l[i] = 16/3 * (gamma_l[i])**2 - 20/3 * (gamma_l[i])**3 + 7/3 * (gamma_l[i])**4
    sigma[i] = alpha[0] + alpha[1] * sigma_s[i] + alpha[2] * sigma_l[i]

# dissociation factor
phi = np.array([1, 1,  1, 2*0.96])


mtac = pd.read_excel('fittedPS/fittedPS_TPM.xlsx', index_col=0)

sse = pd.DataFrame(columns=solutes)

for patient in patientlist:
    
    Cp = data[patient]['cp']
    V = data[patient]['V']
    df_cd = data[patient]['df_cd']
    LS = data[patient]['LS']
    predicted_cd = pd.DataFrame(columns= solutes, dtype = float)
    predicted_V = np.zeros(len(V))
    predicted_cd.loc[0]=df_cd.loc[0]
    predicted_V[0] = V[0]
    x = mtac.loc[patient][:5]
    x[:4]/=(0.998*abya0_s+0.002*abya0_l)
    err, predicted_cd = objective(x, predicted_cd, predicted_V, Cp, V, df_cd, LS)
    sse.loc[patient] = err
    
with pd.ExcelWriter('fittedPS/sse_TPM.xlsx', engine='openpyxl') as writer:
    sse.to_excel(writer)