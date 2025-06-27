# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:09:25 2024

@author: P70073624
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:46:52 2023
Garred model fitting
@author: P70073624
"""



import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import time

MTAC_G = pd.Series(dtype=float) 
solutes = ["Urea", "Creatinine", "Glucose", "Sodium"]
datafile = pd.read_excel('PD_in_silico_csv_export_/PD_in_silico_export_20240723.xlsx')
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)
patientlist = ['PDPC_059', 'PDPC_236','PDPC_268','PDPC_541','PDPC_572','PDPC_800','PDPC_835','PDPC_890','PDPC_946']

mtac = pd.read_excel('fittedPS/fittedPS_GM.xlsx', index_col = 0)

t = 240
sse = pd.DataFrame(columns=solutes)
for i,ind in enumerate(patientlist):
    st = time.time()
    df_cr = datafile.loc[ind,['Creatinine_t20', 'Creatinine_t120', 'Creatinine_t240']]*0.001
    df_cr.index = [20, 120, 240]
    
    df_glu = datafile.loc[ind,['Gluc_t20', 'Gluc_t120', 'Gluc_t240']]
    df_glu.index = [20, 120, 240]
    
    df_ur = datafile.loc[ind,['Ur_t20', 'Ur_t120', 'Ur_t240']]
    df_ur.index = [20, 120, 240]
    
    df_na = datafile.loc[ind,['Na20', 'Na120', 'Na240']]
    df_na.index = [20, 120, 240]
    
    df_cd = pd.concat([df_ur, df_cr, df_glu, df_na], axis = 1)
    df_cd.columns = solutes #mmol/L
    
    
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
    cp = cp*1/(0.984 - 0.000718*70) #70 g/L is mean human total protein blood levels
    
    # Volume from datafiles
    df_V = np.array(datafile.loc[ind,['Vin_t0', 'Vdrain_t240']]) #mL
    df_V[0] += RV0 #assuming a residual volume of 200 ml
    df_V[1] += RV240#assuming a residual volume of 200 ml
    #Linear interpolation to find values of V at all times
    f_V = interp1d([0,240], df_V)
    interpolated_V = f_V(range(0,241))

    V = interpolated_V
    
    x = np.array(mtac.loc[ind][0:4])

    predicted_cd = pd.DataFrame(columns= solutes, dtype = float)
    predicted_cd.loc[0]=df_cd.loc[0]
    for t in range(1,241):
        predicted_cd.loc[t] = cp-(cp-predicted_cd.loc[0])*(V[0]/V[t])*np.exp(-x/V.mean()*t)
        
    df2 = ((df_cd/cp-predicted_cd.loc[df_cd.index]/cp)**2)
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')        
    err = np.sqrt(df2.sum()/len(df_cd))

    sse.loc[ind] = err
    
with pd.ExcelWriter('fittedPS/sse_GM.xlsx', engine='openpyxl') as writer:
    sse.to_excel(writer)  



