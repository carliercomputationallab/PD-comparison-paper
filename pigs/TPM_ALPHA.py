# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:47:16 2024

@author: P70073624
"""


from values import input_values
import multiprocessing
import numpy as np
import random
import scipy
import pandas as pd
import time
import os
import matplotlib.pyplot as plt
from fnmatch import fnmatch



#%%
def objective(x, predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l):
    '''The objective function needed to be minimised'''
    

    t = 240

    predicted_cd, jv = rk(t, x, predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l)

    df2 = ((1-predicted_cd.loc[df_cd.index]/df_cd)**2)
    
    return (sum(np.sqrt(df2.sum()/len(df_cd))), predicted_cd, jv)

#%% OBJECTIVE FUNCTION FOR MINIMISE CANNOT HAVE MORE THAN ONE RETURN OUTPUT, so we 
# repeat the objective function again with just one return
def objective_fn(x, predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l):
    '''The objective function needed to be minimised'''

    t = 240
    
    "call runge Kutta"
    predicted_cd, jv = rk(t, x, predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l)
    nan_mask = np.isnan(predicted_cd)
    predicted_cd = pd.DataFrame(np.where(nan_mask, 0, predicted_cd), columns = solutes)
    #rmse_sodium = ((df_cd['Sodium']/Cp.mean(axis = 0)[2]-predicted_cd.loc[df_cd.index, 'Sodium']/Cp.mean(axis = 0)[2])**2).sum()
    df2 = ((1-predicted_cd.loc[df_cd.index, 'Sodium']/df_cd['Sodium'])**2)
    return np.abs(df2).sum()/len(df_cd)
#%%
#Runge-Kutta
def rk(t, x, predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l):

    Jv = []
    
    for timestep in range(0,t): 
        
        cd = predicted_cd.loc[timestep]
        
        "Apply Runge Kutta Formulas to find next value of y"
        k1, v1 = comdxdt(cd, timestep, x,  predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l)
        k2, v2 = comdxdt(cd + 0.5  *k1, timestep, x,  predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l)
        k3, v3 = comdxdt(cd + 0.5  *k2, timestep, x,  predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l)
        k4, v4 = comdxdt(cd + k3, timestep, x,  predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l)
        
        # Update next value of y
        cd = cd + (1 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        Jv.append( (1 / 6.0)*(v1 + 2 * v2 + 2 * v3 + v4))
        #print(UF)
        predicted_cd.loc[timestep+1] = cd          
    return (predicted_cd, Jv)

# the differential equations
def comdxdt(cd, t, x, predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l):
    
    solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
    
    PS_s[2] = x[0]
    alpha_c = x[1]
    alpha_l = 0.1-alpha_c
    
    #fraction of peritoneal membrane in contact with the dialysis fluid
    af = 16.18 * (1 - np.exp (-0.00077*V[t]))/13.3187
    
    #hydrostatic pressure difference
    delP = delP0 - ((V[t] - (V_fill+Vr))/490)
    
    #peritoneal concentration gradient
    pr = [phi[i]*RT * (Cp[t][i]-cd.iloc[i]) for i in range(len(solutes))]
    sigmas_pr = sum([ sigma_s[i]*pr[i] for i in range(len(solutes))])
    sigmal_pr = sum([ sigma_l[i]*pr[i] for i in range(len(solutes))])

    # #volumetric flows across the pores
    J_vC = af*alpha_c*LS*(delP - sum(pr)) #ml/min
    J_vS = af*alpha[1]*LS*(delP  - sigmas_pr) #ml/min
    J_vL = af*alpha_l*LS*(delP - sigmal_pr) #ml/min

    # #Peclet numbers
    Pe_s = np.array([J_vS  * (1 - sigma_s[i])/(af*PS_s[i]) for i in range(len(solutes))])
    Pe_s[np.isinf(Pe_s)] = 0
    Pe_l = np.array([J_vL  * (1 - sigma_l[i])/(af*PS_l[i]) for i in range(len(solutes))])
    # print(Pe_l, Pe_s)
    
    # #solute flow rate
    J_sS = (J_vS*(1-sigma_s)*(Cp[t]-cd*np.exp(-Pe_s))/(1-np.exp(-Pe_s))).ravel()
    J_sL = (J_vL*(1-sigma_l)*(Cp[t]-cd*np.exp(-Pe_l))/(1-np.exp(-Pe_l))).ravel()

    dxdt = ((J_sS + J_sL)/V[t]-np.array(cd)*(J_vC + J_vS + J_vL-L)/V[t]).ravel()

    return (dxdt, J_vC+J_vS+J_vL)


    
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
files_to_ignore = ['P14', 'P15', 'P18', 'P19', 'P22', 'P23']
patientlist = [
    file for file in patientlist
    if not any(ignore in file for ignore in files_to_ignore)]
# patientlist = ['patient_files/pig2/session1/P10.csv']

# ultrafiltration coefficient
LS = 0.074 #ml/min/mmHg

# fractional pore coefficients
alpha = [0.02, 0.900, 0.08]

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
L = 0.7 #ml/min 

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
phi = np.array([1, 1, 2*0.96, 1, 1, 1])


#%%

if __name__ == '__main__':
    
    ALPHA = pd.DataFrame(columns = ['patient', 'MTAC_sodium' ,'ALPHA-ultrasmall', 'objective fn'], index = range(0, len(patientlist)))
    for i, pfile in enumerate(patientlist):
    
        df_OV = pd.read_excel('fittedPS/fittedPS_TPM.xlsx', index_col = 0).loc[pfile[28:]]
        # from model 7 results for all patients
        MTAC = np.array(df_OV.iloc[:6], dtype=float) / (0.998 * abya0_s + 0.002 * abya0_l)
        PS_s = MTAC * 0.998 * abya0_s 
        PS_l = MTAC * 0.002 * abya0_l
        predicted_cd, Cp, V, df_cd, Vr, V_fill = input_values(pfile)
        V = V*1000
        optimised_values = []
        obj_fn = []
        Jv = np.empty((10,240))
        
        for var in range(5):
            
            #Define initial initial_guess
            x0 = [random.random()* 10, random.random()/10]
            
            '''SLSQP optimisation'''
            result = scipy.optimize.minimize(objective_fn, x0, args = (predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l),
                    method='SLSQP', bounds = [(0, 200), (0, 0.1)],
                    options = {"maxiter" : 1000, "disp": True})
            
            #gather all optimised values
            optimised_values.append(result['x'])
            obj_fn.append(result['fun'])
        ALPHA.loc[i,'patient'] = patientlist[0][28:]
        ALPHA.loc[i, 'MTAC_sodium'] = optimised_values[obj_fn.index(min(obj_fn))][0]
        ALPHA.loc[i, 'ALPHA-ultrasmall'] = optimised_values[obj_fn.index(min(obj_fn))][1]
        ALPHA.loc[i,'objective fn'] = min(obj_fn)
        
        predicted_cd, Cp, V, df_cd, Vr, V_fill = input_values(pfile)
        
        V = V*1000
        
        _ , predicted_cd, _ = objective(optimised_values[obj_fn.index(min(obj_fn))], predicted_cd, Cp, V, df_cd, Vr, V_fill, PS_s, PS_l)
        fig, ax = plt.subplots(3,2, figsize = (12,18))
        t = 240
        #urea
        df_cd['Urea'].plot( ax = ax[0,0], label = 'data', style = '.')
        ax[0,0].plot(np.arange(t+1),predicted_cd['Urea'], label = 'predicted')
        # ax[0,0].text(0.6, 0.1, f'MTAC = {result["x"][0]:.2f} ml/min', transform=ax[0,0].transAxes)
        ax[0,0].set_title("Urea")
        
        #creatinine
        df_cd['Creatinine'].plot( ax = ax[0,1], style = '.')
        ax[0,1].plot(np.arange(t+1),predicted_cd['Creatinine'])
        # ax[0,1].text(0.6, 0.1, f'MTAC = {result["x"][1]:.2f} ml/min', transform=ax[0,1].transAxes)
        ax[0,1].set_title("Creatinine")
        
        #Sodium
        df_cd['Sodium'].plot( ax = ax[1,0],  style = '.')
        ax[1,0].plot(np.arange(t+1),predicted_cd['Sodium'] )
        # ax[1,0].text(0.6, 0.5, f'MTAC = {result["x"][2]:.2f} ml/min', transform=ax[1,0].transAxes)
        ax[1,0].set_title("Sodium")
        
        #Phosphate
        df_cd['Phosphate'].plot( ax = ax[1,1], style = '.')
        ax[1,1].plot(np.arange(t+1),predicted_cd['Phosphate'] )
        # ax[1,1].text(0.6, 0.1, f'MTAC = {result["x"][3]:.2f} ml/min', transform=ax[1,1].transAxes)
        ax[1,1].set_title("Phosphate")
        
        #Glucose
        df_cd['Glucose'].plot( ax = ax[2,0], style = '.')
        ax[2,0].plot(np.arange(t+1),predicted_cd['Glucose'])
        # ax[2,0].text(0.6, 0.5, f'MTAC = {result["x"][4]:.4f} ml/min', transform=ax[2,0].transAxes)
        ax[2,0].set_title("Glucose")
        
        #Potassium
        df_cd['Potassium'].plot( ax = ax[2,1], style = '.')
        ax[2,1].plot(np.arange(t+1),predicted_cd['Potassium'])
        # ax[2,1].text(0.6, 0.1, f'MTAC = {result["x"][5]:.2f} ml/min', transform=ax[2,1].transAxes)
        ax[2,1].set_title("Potassium")
        
        fig.supxlabel("time, min")
        fig.supylabel("Dialysate concentration, mmol")
        plt.title(f'{pfile[28:]}')
        plt.subplots_adjust(top=0.88,
                            bottom=0.11,
                            left=0.09,
                            right=0.9,
                            hspace=0.295,
                            wspace=0.215)
        with pd.ExcelWriter(f'fittedPS/corrected_predicted_cd_TPM_{pfile[28:-4]}.xlsx', engine='openpyxl') as writer:
            predicted_cd.to_excel(writer) 
    with pd.ExcelWriter('fittedPS/fittedPS_alpha_TPM.xlsx', engine='openpyxl') as writer:
        ALPHA.to_excel(writer) 
    #%%
    
    # "Use this to plot"
    
    
    # predicted_cd, Cp, V, df_cd, Vr, V_fill = input_values(patientlist[0])
    
    # V = V*1000
    
    # _ , predicted_cd, _ = objective(optimised_values[obj_fn.index(min(obj_fn))], predicted_cd, Cp, V, df_cd, Vr, V_fill)
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
    # plt.subplots_adjust(top=0.88,
    #                     bottom=0.11,
    #                     left=0.09,
    #                     right=0.9,
    #                     hspace=0.295,
    #                     wspace=0.215)
    
    # # with pd.ExcelWriter('predicted_files/corrected_predicted_cd_TPM_p10.xlsx', engine='openpyxl') as writer:
    # #     predicted_cd.to_excel(writer) 
