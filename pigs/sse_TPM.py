# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 12:35:35 2023

@author: P70073624
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 11:35:19 2022

@author: P70073624
"""


from values import input_values
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from fnmatch import fnmatch



#%%
def objective(x, predicted_cd, predicted_V, Cp, V, df_cd, Vr, V_fill):
    '''The objective function needed to be minimised'''
    

    t = 240

    predicted_cd,_, jv = rk(t, x, predicted_cd, predicted_V, Cp, df_cd, Vr, V_fill)

    
    df2 = ((df_cd/Cp[0]-predicted_cd.loc[df_cd.index]/Cp[0])**2).astype('float64')
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
    return np.sqrt(df2.sum()/len(df_cd)), predicted_cd



#%%
#Runge-Kutta
def rk(t, x, predicted_cd, predicted_V, Cp, df_cd, Vr, V_fill):

    Jv = []
    
    for timestep in range(0,t): 
        
        cd = np.array(predicted_cd.loc[timestep])
        Vp = predicted_V[timestep]
        V0 = predicted_V[0]
        "Apply Runge Kutta Formulas to find next value of y"
        k1, v1 = comdxdt(cd, timestep, x,  predicted_cd, Cp, Vp, df_cd, V0)
        k2, v2 = comdxdt(cd + 0.5  *k1, timestep, x,  predicted_cd, Cp, Vp + v1/2, df_cd, V0)
        k3, v3 = comdxdt(cd + 0.5  *k2, timestep, x,  predicted_cd, Cp, Vp + v2/2, df_cd, V0)
        k4, v4 = comdxdt(cd + k3, timestep, x,  predicted_cd, Cp, Vp+v3, df_cd, V0)
        
        # Update next value of y
        cd = cd + (1 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        Jv.append( (1 / 6.0)*(v1 + 2 * v2 + 2 * v3 + v4))
        #print(UF)
        predicted_cd.loc[timestep+1] = cd     
        predicted_V[timestep+1] = predicted_V[timestep] + (1 / 6.0)*(v1 + 2 * v2 + 2 * v3 + v4) 
    return (predicted_cd, predicted_V, Jv)

# the differential equations
def comdxdt(cd, t, x, predicted_cd, Cp, Vp, df_cd, V0):
    
    solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
    
    #MTAC for small pores
    PS_s = x[0:6] * 0.998 * abya0_s #fraction small pore surface area - Rippe, A THREE-PORE MODEL OF PERITONEAL TRANSPORT table 1

    #MTAC for large pores
    PS_l = x[0:6] * 0.002 * abya0_l# Ref: two pore model-oberg,rippe, table 1, A0L/A0 value
    
    # alpha[0] = x[6]
    # alpha[1] = 0.92-x[6]
    # alpha[2] = 0.08
    L = x[6]
    
    #fraction of peritoneal membrane in contact with the dialysis fluid
    af = 16.18 * (1 - np.exp (-0.00077*Vp))/13.3187
    
    #hydrostatic pressure difference
    delP = delP0 - ((V[t] - V0)/490)
    
    #peritoneal concentration gradient
    pr = [phi[i]*RT * (Cp[t][i]-cd[i]) for i in range(len(solutes))]
    sigmas_pr = sum([phi[i] * sigma_s[i]*pr[i] for i in range(len(solutes))])
    sigmal_pr = sum([phi[i] * sigma_l[i]*pr[i] for i in range(len(solutes))])

    # #volumetric flows across the pores
    J_vC = af*alpha[0]*LS*(delP - sum(pr)-22) #ml/min
    J_vS = af*alpha[1]*LS*(delP  - sigmas_pr-0.963*22) #ml/min
    J_vL = af*alpha[2]*LS*(delP - sigmal_pr-0.083*22) #ml/min

    # #Peclet numbers
    Pe_s = np.array([J_vS  * (1 - sigma_s[i])/(af*PS_s[i]) for i in range(len(solutes))])
    Pe_l = np.array([J_vL  * (1 - sigma_l[i])/(af*PS_l[i]) for i in range(len(solutes))])
    
    # #solute flow rate
    J_sS = (J_vS*(1-sigma_s)*(Cp[t]-cd*np.exp(-Pe_s))/(1-np.exp(-Pe_s))).ravel()
    J_sL = (J_vL*(1-sigma_l)*(Cp[t]-cd*np.exp(-Pe_l))/(1-np.exp(-Pe_l))).ravel()

    dxdt = ((J_sS + J_sL)/Vp-np.array(cd)*(J_vC + J_vS + J_vL-L)/Vp).ravel()

    return (dxdt, J_vC+J_vS+J_vL)

    

# ultrafiltration coefficient
LS = 0.074 #ml/min/mmHg

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

solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
#radius of molecules
r = np.array([ 2.6, 3.0, 2.3, 2.77, 3.7, 2.8]) #the phosphate radius is approximated from its topological surface area
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
phi = np.array([1, 1, 2*0.96, 1, 1, 1])
#%%

if __name__ == '__main__':
    
    
    #%%
    # for a 7:4 training to test split
    "Get all MTACs for the pig in question"
    root = 'patient_files/'
    pattern = "*.csv"
    patientlist = []
    for path, subdirs, files in os.walk(root):
        for name in files:
            if fnmatch(name, pattern):
                patientlist.append(os.path.join(path, name))#random.sample(os.listdir('patient_files/pig2/session1'),7) #randomly generated using random.sample

    # patientlist = ['patient_files/pig2/session1/P10.csv']
    mtac = pd.read_excel('fittedPS/fittedPS_TPM.xlsx', index_col=0)
    
    sse = pd.DataFrame(columns=solutes)
    
    for patient in patientlist:
        predicted_cd, Cp, V, df_cd, Vr, V_fill, _ = input_values(patient)
        predicted_V = np.zeros(len(V))
        predicted_V[0] = V[0]
        patient = patient[28:]
        x = mtac.loc[patient][:7]
        print(x)
        x[0:6] = x[0:6]/(0.998*abya0_s+0.002*abya0_l)
        err, predicted_cd = objective(x, predicted_cd, predicted_V, Cp, V, df_cd, Vr, V_fill)
        sse.loc[patient] = err
        
    with pd.ExcelWriter('fittedPS/sse_TPM.xlsx', engine='openpyxl') as writer:
        sse.to_excel(writer)
    
    # with pd.ExcelWriter('fittedPS/predicted_cd_m7.xlsx', engine='openpyxl') as writer:
    #     predicted_cd.to_excel(writer)
    
    # # "Use for plot"
    # fig, ax = plt.subplots(3,2, figsize = (12,18))
    # t = 240
    # #urea
    # df_cd['Urea'].plot( ax = ax[0,0], label = 'data', style = '.')
    # ax[0,0].plot(np.arange(t+1),predicted_cd['Urea'], label = 'predicted')
    # # ax[0,0].text(0.6, 0.1, f'MTAC = {result["x"][0]:.2f} ml/min', transform=ax[0,0].transAxes)
    # ax[0,0].set_title("Urea")

    # #creatinine
    # df_cd['Creatinine'].plot( ax = ax[0,1], style = '.')
    # ax[0,1].plot(np.arange(t+1),predicted_cd['Creatinine'])
    # # ax[0,1].text(0.6, 0.1, f'MTAC = {result["x"][1]:.2f} ml/min', transform=ax[0,1].transAxes)
    # ax[0,1].set_title("Creatinine")

    # #Sodium
    # df_cd['Sodium'].plot( ax = ax[1,0],  style = '.')
    # ax[1,0].plot(np.arange(t+1),predicted_cd['Sodium'] )
    # # ax[1,0].text(0.6, 0.5, f'MTAC = {result["x"][2]:.2f} ml/min', transform=ax[1,0].transAxes)
    # ax[1,0].set_title("Sodium")

    # #Phosphate
    # df_cd['Phosphate'].plot( ax = ax[1,1], style = '.')
    # ax[1,1].plot(np.arange(t+1),predicted_cd['Phosphate'] )
    # # ax[1,1].text(0.6, 0.1, f'MTAC = {result["x"][3]:.2f} ml/min', transform=ax[1,1].transAxes)
    # ax[1,1].set_title("Phosphate")

    # #Glucose
    # df_cd['Glucose'].plot( ax = ax[2,0], style = '.')
    # ax[2,0].plot(np.arange(t+1),predicted_cd['Glucose'])
    # # ax[2,0].text(0.6, 0.5, f'MTAC = {result["x"][4]:.4f} ml/min', transform=ax[2,0].transAxes)
    # ax[2,0].set_title("Glucose")

    # #Potassium
    # df_cd['Potassium'].plot( ax = ax[2,1], style = '.')
    # ax[2,1].plot(np.arange(t+1),predicted_cd['Potassium'])
    # # ax[2,1].text(0.6, 0.1, f'MTAC = {result["x"][5]:.2f} ml/min', transform=ax[2,1].transAxes)
    # ax[2,1].set_title("Potassium")

    # fig.supxlabel("time, min")
    # fig.supylabel("Dialysate concentration, mmol")
    # plt.suptitle("Model 1 predictions of dialysate concentration")
    # plt.subplots_adjust(top=0.88,
    #                     bottom=0.11,
    #                     left=0.09,
    #                     right=0.9,
    #                     hspace=0.295,
    #                     wspace=0.215)