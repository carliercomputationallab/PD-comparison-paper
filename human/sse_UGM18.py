# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:52:26 2024

@author: P70073624
"""




import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



#%%
def calculate_penalty_fct(x):
    penalty = [0, 0, 0, 0]
    for param in x[4:8]:
        if param < 0:
            penalty[i]+=param ** 2
            
        elif param > 1:
            penalty[i]+=((param-1)**2)       
    return penalty
            
def calculate_penalty_sico(x):
    penalty = [0, 0, 0, 0]
    for i,param in enumerate(x[8:12]):
        if param < 0:
            penalty[i]+=param ** 2
            
        elif param > 1:
            penalty[i]+=((param-1)**2)        
    return penalty       

def objective(x, predicted_cd, Cp, V, df_cd):
    '''The objective function needed to be minimised'''

    t = 240
    predicted_cd = rk(t, x, predicted_cd, Cp, V, df_cd)
    df2 = ((df_cd/Cp-predicted_cd.loc[df_cd.index]/Cp)**2)
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
    cp = Cp # plasma solute concentration at time t
    #see models of graff for explanation  
    MTAC = np.array(x[0:4])
    fct = np.array(x[4:8])
    fct = np.around(fct, decimals=4)
    SiCo = np.array(x[8:12]) #SiCo
    SiCo = np.around(SiCo, decimals=4)
    L = 0.3 # lymphatic flow rate
    QU = V[t+1] - V[t] + L #UF rate
    beta = QU * SiCo/MTAC.ravel() #ml/min -> l/min
    f = np.array([0 if b > 30 else 0.5 if b == 0 else (1/b)-1/(np.exp(b)-1) for b in beta]) #check model 2 for explanation
    MC = cp-f*(cp-cd)
    Cl = cp if L < 0 else cd
    dxdt = (MTAC * (fct * cp - cd) + SiCo * QU * MC - L * Cl - cd * QU)/V[t].ravel()
    
    
    
    return dxdt


    
datafile = pd.read_excel('PD_in_silico_csv_export_/PD_in_silico_export_20240723.xlsx')
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)
patientlist = ['PDPC_059', 'PDPC_236','PDPC_268','PDPC_541','PDPC_572','PDPC_800','PDPC_835','PDPC_890','PDPC_946']


#create a dictionary to save all patient data
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
    RV0 = datafile.loc[ind, "RV0"]
    RV240 = datafile.loc[ind, "RV240"]
        
    print(RV0, RV240, flush= True)
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
    print(df_V, flush = True)
    df_V[0] += RV0 #assuming a residual volume of 200 ml
    df_V[1] += RV240
    #Linear interpolation to find values of V at all times
    f_V = interp1d([0,240], df_V)
    interpolated_V = f_V(range(0,241))

    V = interpolated_V

    # Set the ultrafiltraiton coefficients
    LS = 0.074 #ml/min mmHg

    pat_d = {'df_cd': df_cd, 'cp': cp, 'V': V,  'LS':LS}
    
    data[ind] = pat_d
  


mtac = pd.read_excel('fittedPS/fittedPS_UGM18.xlsx', index_col = 0)

solutes = ["Urea", "Creatinine",  "Glucose", "Sodium"]

sse = pd.DataFrame(columns=solutes)

for patient in patientlist:
    Cp = data[patient]['cp']
    V = data[patient]['V']
    df_cd = data[patient]['df_cd']
    predicted_cd = pd.DataFrame(columns= ["Urea", "Creatinine", "Glucose", "Sodium"], dtype = float)
    predicted_cd.loc[0]=df_cd.loc[0]
    x = mtac.loc[patient][:12]
    err, predicted_cd = objective(x, predicted_cd, Cp, V, df_cd)
    sse.loc[patient] = err
    
with pd.ExcelWriter('fittedPS/sse_UGM18.xlsx', engine='openpyxl') as writer:
    sse.to_excel(writer)
    
    
