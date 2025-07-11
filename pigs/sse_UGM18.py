# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 15:50:07 2023

@author: P70073624
"""


from values import input_values
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from fnmatch import fnmatch



#%%
def calculate_penalty_fct(x):
    penalty = [0] * 6
    for i, param in enumerate(x[6:12]):
        if param < 0:
            penalty[i] = param ** 2
            
        elif param > 1:
            penalty[i] = (param-1)**2
            
    return penalty
            
def calculate_penalty_sico(x):
    penalty = [0] * 6
    for i, param in enumerate(x[12:18]):
        if param < 0:
            penalty[i] = param ** 2
            
        elif param > 1:
            penalty[i] = (param-1)**2
            
    return penalty        

def objective(x, predicted_cd, Cp, V, df_cd):
    '''The objective function needed to be minimised'''

    t = 240
    predicted_cd = rk(t, x, predicted_cd, Cp, V, df_cd)
    df2 = ((df_cd/Cp[0]-predicted_cd.loc[df_cd.index]/Cp[0])**2).astype('float64')
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
    err = np.sqrt(df2.sum()/len(df_cd))
    penalty = np.add(calculate_penalty_fct(x),calculate_penalty_sico(x))
    return [sum(x) for x in zip(err, penalty)], predicted_cd


#%%
#Runge-Kutta
def rk(t, x, predicted_cd, Cp, V, df_cd):
    
        
    for timestep in range(0,t): 
        
        cd = predicted_cd.loc[timestep]
        
        "Apply Runge Kutta Formulas to find next value of y"
        k1 = compute(cd, timestep, x,  predicted_cd, Cp, V, df_cd)
        k2 = compute(cd + 0.5  *k1, timestep, x,  predicted_cd, Cp, V, df_cd)
        k3 = compute(cd + 0.5  *k2, timestep, x,  predicted_cd, Cp, V, df_cd)
        k4 = compute(cd + k3, timestep, x,  predicted_cd, Cp, V, df_cd)

        # Update next value of y
        cd = cd + (1 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)

        predicted_cd.loc[timestep+1] = cd
        
    return predicted_cd


# the differential equations
def compute(cd, t, x, predicted_cd, Cp, V, df_cd):
    '''
    

    Parameters
    ----------
    cd : predicted dialysate concentration
    t : timepoint
        DESCRIPTION.
    x : intial matrix
        x[0:6] = MTAC
        x[6:12] = fct
        x[12:18] = SiCo
        x[18] = QL
    model : 1-6
        DESCRIPTION.

    Returns
    -------
    dxdt : conc gradient
        derivative

    '''
    
    # solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
    cp = Cp[t] # plasma solute concentration at time t
    #see models of graff for explanation  
    MTAC = np.array(x[0:6])
    fct = np.array(x[6:12])
    SiCo = np.array(x[12:18]) #SiCo
    L = 0.3 # lymphatic flow rate
    QU = V[t+1] - V[t] + L #UF rate
    beta = QU * SiCo/MTAC.ravel() #ml/min -> l/min
    f = np.array([0 if b > 30 else 0.5 if b == 0 else (1/b)-1/(np.exp(b)-1) for b in beta]) #check model 2 for explanation
    MC = cp-f*(cp-cd)
    Cl = cp if L < 0 else cd
    dxdt = (MTAC * (fct * cp - cd) + SiCo * QU * MC - L * Cl - cd * QU)/V[t].ravel()
    
    
    return dxdt


    
# patientlist= random.sample(os.listdir('patient_files/pig2/session1'),7) #randomly generated using random.sample
"Get all MTACs for the pig in question"
root = 'patient_files/'
pattern = "*.csv"
patientlist = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            patientlist.append(os.path.join(path, name))
# patientlist = [p for p in patientlist if 'pig2' not in p] #since we already have the results from pig 2
# patientlist = ['patient_files/pig2/session1/P10.csv']

    
    # df_OV.loc['mean'] = df_OV[df_OV['objective function']<1e5].mean()
mtac = pd.read_excel('fittedPS/fittedPS_UGM18.xlsx', index_col = 0)

solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]

sse = pd.DataFrame(columns=solutes)

for patient in patientlist:
    predicted_cd, Cp, V, df_cd, Vr, V_fill, _ = input_values(patient)
    patient = patient[28:]
    x = mtac.loc[patient][:18]
    err, predicted_cd = objective(x, predicted_cd, Cp, V, df_cd)
    # "Use for plot"
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
    sse.loc[patient] = err
    
with pd.ExcelWriter('fittedPS/sse_UGM18.xlsx', engine='openpyxl') as writer:
    sse.to_excel(writer)
