# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 12:08:54 2024

@author: P70073624
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:46:52 2023
simple Waniewski model fitting
@author: P70073624
"""



import time
import numpy as np
import os
import pandas as pd
import random
from fnmatch import fnmatch
from scipy.interpolate import interp1d
from sklearn.preprocessing import MinMaxScaler


MTAC_W = pd.Series(dtype=float) 
solutes = ["Urea", "Creatinine", "Glucose"]

#get patient data
datafile = pd.read_excel('PD_in_silico_csv_export_/PD_in_silico_export_20240723.xlsx')
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)
patientlist = ['PDPC_059', 'PDPC_236','PDPC_268','PDPC_541','PDPC_572','PDPC_800','PDPC_835','PDPC_890','PDPC_946']

comp_time = []
obj_fn = []
t = 240
for i,pfile in enumerate(patientlist):
    st = time.time()
    
    # dialysate concentrations from data file
    df_cr = datafile.loc[pfile,['Creatinine_t20', 'Creatinine_t120', 'Creatinine_t240']]*0.001
    df_cr.index = [20, 120, 240]
    
    df_glu = datafile.loc[pfile,['Gluc_t20', 'Gluc_t120', 'Gluc_t240']]
    df_glu.index = [20, 120, 240]
    
    df_ur = datafile.loc[pfile,['Ur_t20', 'Ur_t120', 'Ur_t240']]
    df_ur.index = [20, 120, 240]
    
    df_na = datafile.loc[pfile,['Na20', 'Na120', 'Na240']]
    df_na.index = [20, 120, 240]
    
    df_cd = pd.concat([df_ur, df_cr, df_glu, df_na], axis = 1)
    df_cd.columns = ["Urea", "Creatinine", "Glucose", "Sodium"] #mmol/L
    
    #fill volume from datafile    
    Vfill = datafile.loc[pfile, 'Vin_t0']

    
    RV0 = datafile.loc[pfile, "RV0"]
    RV240 = datafile.loc[pfile, "RV240"]
        
    print(RV0, RV240, flush= True)
    ur_0 = [datafile.loc[pfile, 'Ur_night']*RV0/(Vfill+RV0) if datafile.loc[pfile, 'Ur_night'] > 0 else 0]
    cr_0 = [datafile.loc[pfile, 'Cr_night']*RV0/(Vfill+RV0)*0.001 if datafile.loc[pfile, 'Cr_night'] > 0 else 0]
    glu_0 = [214 if datafile.loc[pfile, 'Gluc_PET'] == 3.86 else  126 if datafile.loc[pfile, 'Gluc_PET'] == 2.27 else  75 ]
    na_0 = [132]
    df_cd.loc[0] = np.concatenate((ur_0, cr_0, glu_0, na_0))
    

    #blood plasma concentration from datafiles
    cp = datafile.loc[pfile, ['Ur_t2_B','Cr_t2_B', 'Gluc_t2_B', 'Na_t2_B']]
    cp.index = ["Urea", "Creatinine", "Glucose", "Sodium"]
    if cp.loc['Glucose'] < 0:
        cp.loc['Glucose'] = 6.0
        
    cp.loc['Creatinine'] *=0.001
    cp.loc['Sodium'] *=0.94
    cp = cp*1/(0.984 - 0.000718*70) #70 g/L is mean human total protein blood levels
    
    # Volume from datafiles
    df_V = np.array(datafile.loc[pfile,['Vin_t0', 'Vdrain_t240']]) #mL
    df_V[0] += RV0 #assuming a residual volume of 200 ml
    df_V[1] += RV240
    f_V = interp1d([0,240], df_V)
    interpolated_V = f_V(range(0,241))
    V = interpolated_V
    
    #calculate MTAC based on Waniewski model
    mtac = V.mean()/t*np.log(
        pd.to_numeric((pow(V[0],(1/2))*(cp-df_cd.loc[0]))/
        (pow(V[240],(1/2))*(cp-df_cd.loc[240]))))
    
    # if the mtac is NAN by some chance, try to calculate MTAC again by using another data point
    for solute in solutes:
        for ind in list(df_cd.index)[::-1]:
            mtac[solute] = V.mean()/ind*np.log(
                (pow(V[0], 0.5)*(cp[solute]-df_cd[solute].loc[0]))/
                (pow(V[ind], 0.5)*(cp[solute]-df_cd[solute].loc[ind])))
            if np.isfinite(mtac[solute]):
                break
    
    MTAC_W = pd.concat([MTAC_W,mtac], axis=1)
    
    
    predicted_cd = pd.DataFrame(columns= ["Urea", "Creatinine", "Glucose", "Sodium"], dtype = float)
    predicted_cd.loc[0]=df_cd.loc[0]
    
    
    # calculate predicted dialysate cocnetration
    for t in range(1,241):
        predicted_cd.loc[t] = cp-(cp-predicted_cd.loc[0])*((V[0]/V[t])**0.5)*np.exp(-mtac/V.mean()*(t))
    predicted_cd.to_excel(f'predicted_cd/SWM_{pfile}.xlsx')     
    #error calculation
    df2 = ((df_cd/cp-predicted_cd.loc[df_cd.index]/cp)**2)
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
    
    
    #calculate RMSE
    obj_fn.append(sum(np.sqrt(df2.sum()/len(df_cd))))            
    et = time.time()
    comp_time.append(et-st)
    # with pd.ExcelWriter('fittedPS/predicted_cd_SWM_PDPC_754.xlsx', engine='openpyxl') as writer:
    #     predicted_cd.to_excel(writer) 
    
MTAC_W = MTAC_W.T
MTAC_W = MTAC_W[1:]
MTAC_W.index = patientlist
MTAC_W.columns = ['MTAC_urea', 'MTAC_crea', 'MTAC_glu', 'MTAC_sodium']
MTAC_W['obj_fn'] = obj_fn
MTAC_W['comp-time'] = comp_time
print(MTAC_W)



with pd.ExcelWriter('fittedPS/fittedPS_SWM.xlsx', engine='openpyxl') as writer:
    MTAC_W.to_excel(writer)    



        