

# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 16:23:23 2023
UNIFIED GRAFF MODEL
Check corresponding file in human for details
the only difference is 6 solutes instead of 3.
and the data from the patient files is collected by the file "values"
make sure that is in the same folder.
@author: P70073624
"""

from values import input_values
import multiprocessing
import numpy as np
import scipy
import pandas as pd
import time
import os
import matplotlib.pyplot as plt
from fnmatch import fnmatch
from scipy.optimize import Bounds


#%%
def calculate_penalty(x):
    penalty = 0
    for param in x[6:18]:
        if param < 0:
            penalty += param ** 2
            
        elif param > 1:
            penalty += (param-1)**2
            
    return penalty
            
        

def objective(x, predicted_cd, Cp, V, df_cd):
    '''The objective function needed to be minimised'''

    t = 240
    predicted_cd = rk(t, x, predicted_cd, Cp, V, df_cd)
    df2 = ((df_cd/Cp[0]-predicted_cd.loc[df_cd.index]/Cp[0])**2).astype('float64')
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
    penalty = calculate_penalty(x)
    return (sum(np.sqrt(df2.sum()/len(df_cd))) + penalty, predicted_cd)

#%%
def objective_fn(x, predicted_cd, Cp, V, df_cd):
    '''The objective function needed to be minimised'''
    # print(x)

    t = 240
    predicted_cd = rk(t, x, predicted_cd, Cp, V, df_cd)
    df2 = ((df_cd/Cp[0]-predicted_cd.loc[df_cd.index]/Cp[0])**2).astype('float64')
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
    penalty = calculate_penalty(x)
    
    # predicted_cd[predicted_cd>200] = np.nan
    # print(predicted_cd.loc[df_cd.index])
    try:
        return sum(np.sqrt(df2.sum()/len(df_cd))) + penalty
    except ValueError:
        return 1000.0
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
    fct = np.around(fct, decimals=4)
    SiCo = np.array(x[12:18]) #SiCo
    SiCo = np.around(SiCo, decimals=4)
    L = 0.3 # lymphatic flow rate
    QU = V[t+1] - V[t] + L #UF rate
    beta = QU * SiCo/MTAC.ravel() #ml/min -> l/min
    f = np.array([0 if b > 30 else 0.5 if b == 0 else (1/b)-1/(np.exp(b)-1) for b in beta]) #check model 2 for explanation
    MC = cp-f*(cp-cd)
    Cl = cp if L < 0 else cd
    dxdt = (MTAC * (fct * cp - cd) + SiCo * QU * MC - L * Cl - cd * QU)/V[t].ravel()
    
    
    return dxdt

def multiprocessing_func(pfile):
    
    st = time.time()
    Nx = 18 #6 MTAC, 4 fct (except crea and glucose), 5 SiCo (except glucose), 1 L
    predicted_cd, Cp, V, df_cd, _, _, _ = input_values(pfile)    
    optimised_values = np.empty(Nx)
    obj_fn = []
    print(pfile)
    
    for var in range(10):
        
        #Define initial initial_guess
        x0 = np.concatenate((np.random.randint(1, 21, size=6), np.random.randint(2, size=12)))
        # lbound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5]
        # ubound = [200, 200, 200, 200, 200, 200,  10, 10, 10, 10, 10,  10, 10, 10, 10, 10, 10, 10] 
        # bounds = Bounds(lbound, ubound)
        lbound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ubound = [200, 200, 5, 200, 200, 200,  1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1] 
        bounds = Bounds(lbound, ubound)
        
        '''SLSQP optimisation'''
        result = scipy.optimize.minimize(objective_fn, x0, args = (predicted_cd, Cp, V, df_cd),
                method='SLSQP', bounds = bounds,
                options = {"maxiter" : 1000, "disp": True})
        
        #gather all optimised values
        optimised_values = np.vstack((optimised_values,result['x'].tolist()))
        obj_fn.append(result['fun'])
       
    et = time.time()    
    return (optimised_values, obj_fn, et-st)
    
# patientlist= random.sample(os.listdir('patient_files/pig2/session1'),7) #randomly generated using random.sample
"Get all MTACs for the pig in question"
root = 'patient_files/'
pattern = "*.csv"
patientlist = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            patientlist.append(os.path.join(path, name))

if __name__ == '__main__':
    
    
        
    "parallel processing"
    pool = multiprocessing.Pool(8)
    result_list = pool.map(multiprocessing_func,patientlist)
    pool.close()
    
    
    "single processor"
    # Initialize the Min-Max scaler
    # scaler = MinMaxScaler()
    # result_list = []
    # for pfile in patientlist:
    #     Nx = 15 #6 MTAC, 4 fct (except crea and glucose), 5 SiCo (except glucose), 1 L
    #     predicted_cd, Cp, V, df_cd, _, _ = input_values(pfile)    
    #     optimised_values = np.empty(Nx)
    #     obj_fn = []
    #     print(pfile)
        
    #     for var in range(10):
            
    #         #Define initial initial_guess
    #         x0 = np.concatenate((np.random.randint(1, 21, size=6), np.random.randint(2, size=8), np.array([np.random.uniform(0, 0.00001)])))
    #         lbound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, -5, -5, -5, -5, -5, -5, -5, -5,  -200]
    #         ubound = [200, 200, 200, 200, 200, 200, 10, 10, 10, 10, 10, 10, 10, 10, 200] 
    #         bounds = Bounds(lbound, ubound)
            
    #         '''SLSQP optimisation'''
    #         result = scipy.optimize.minimize(objective_fn, x0, args = (predicted_cd, Cp, V, df_cd),
    #                 method='SLSQP', bounds = bounds,
    #                 options = {"maxiter" : 1000, "disp": False})
            
    #         #gather all optimised values
    #         optimised_values = np.vstack((optimised_values,result['x'].tolist()))
    #         obj_fn.append(result['fun'])
            
    #     result_list.append([optimised_values, obj_fn])
                
                
    
    
   
    
    
    #%%
    OF = []
    OV = np.empty(len(patientlist),dtype=object)
    time_TPM = np.empty(len(patientlist),dtype=object)
    
    for i in range(len(result_list)):
        OF.append(min(result_list[i][1]))
        OV[i] = result_list[i][0][np.argmin(result_list[i][1])+1]
        time_TPM[i] = result_list[i][2] #collect all volume flow through pores
        predicted_cd, Cp, V, df_cd, Vr, V_fill, _ = input_values(patientlist[i])
        _, predicted_cd = objective(OV[i], predicted_cd, Cp, V, df_cd) #collect all volume flow through pores
        predicted_cd.to_excel(f'predicted_cd/UGM18_{patientlist[i][28:]}.xlsx') 
        
        
    #%% 
    cols = ['MTAC_urea', 'MTAC_crea','MTAC_sodium', 'MTAC_phosphate','MTAC_glu', 'MTAC_potassium',
            'fct_urea', 'fct_crea','fct_sodium', 'fct_phosphate','fct_glu', 'fct_potassium',
            'sico_urea', 'sico_crea','sico_sodium', 'sico_phosphate', 'sico_glu', 'sico_potassium'] 
    df_OV = pd.DataFrame([arr for arr in OV], columns=cols)   
    df_OV['obj_fn'] = OF
    patient = [lst[28:] for lst in patientlist]
    df_OV.index = patient
    df_OV['comp time'] = time_TPM
    df_OV.to_excel('fittedPS/fittedPS_UGM18.xlsx')


    
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
    