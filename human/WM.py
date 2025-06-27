# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 12:15:44 2024
new waniewski model
@author: P70073624
"""



import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import time
from sklearn.linear_model import LinearRegression
from scipy.integrate import quad
import matplotlib.pyplot as plt

def integrand1(t, df):
    t = int(t)   
    return cp[j] - df.loc[t, solute]

def integrand2(t, df):
    t = int(t)
    return Qv* ((1- F)* cp[j]+ F* df.loc[t, solute])
    
solutes = ["Urea", "Creatinine", "Glucose", "Sodium"]

#%%
#get patient data
datafile = pd.read_excel('PD_in_silico_csv_export_/PD_in_silico_export_20240723.xlsx')
datafile.replace(-99, np.nan, inplace = True)
datafile.set_index('Participant Id', inplace = True)
patientlist = ['PDPC_059', 'PDPC_236','PDPC_268','PDPC_541','PDPC_572','PDPC_800','PDPC_835','PDPC_890','PDPC_946']

comp_time = []
obj_fn = []
data = []
t = 240
for i,ind in enumerate(patientlist):
    st = time.time()
    
    # dialysate concentrations from data file
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

    
    df_cd.sort_index(inplace = True)
    
    index = pd.Series([20,120, 240])
    df_cd = df_cd.set_index([index])
    df_cd = df_cd.astype("float64")
    
    
    #fill volume from datafile    
    Vfill = datafile.loc[ind, 'Vin_t0']
    
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
    
    new_index = pd.Index(np.arange(0, 241, 1))
    df_cd2 = df_cd.copy(deep = True)
    # Reindex the DataFrame
    df_cd2 = df_cd2.reindex(new_index)
    df_cd2 = df_cd2.interpolate(method = 'index', limit_direction = "both")
    

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
    df_V[1] += RV240
    f_V = interp1d([0,240], df_V)
    interpolated_V = f_V(range(0,241))
    V = interpolated_V
#%%    
    # print(df_V)
    
    # define F. F can vary between 0 to 1
    F = 0.5 #dimensionless
    timesteps = np.sort( np.array(df_cd.index))
    X1 = np.empty(len(timesteps)-1)
    X2 = np.empty(len(timesteps)-1)
    VdCd = np.empty(len(timesteps)-1)
    K = np.empty(len(solutes)) #represents MTAC
    S = np.empty(len(solutes)) #represents SiCo
    predicted_cd = pd.DataFrame(columns= solutes, dtype = float)
    
    predicted_cd.loc[0]=df_cd.loc[0] #the value of predicted cd0 comes from the initial guess

    # implementation of Waniewski model
    for j, solute in enumerate(solutes):
        for i in range(1,len(timesteps)):
            X1[i-1] = quad(integrand1, timesteps[i-1], timesteps[i], df_cd2)[0]
            Qv = (V[240]-V[0])/240 - 0.3 *240
            X2[i-1] = quad(integrand2, timesteps[i-1], timesteps[i], df_cd2)[0]
            VdCd[i-1] = V[timesteps[i]]*df_cd.loc[timesteps[i], solute] - V[timesteps[i-1]]*df_cd.loc[timesteps[i-1], solute]

        X1 = X1.reshape(-1,1)
        X2 = X2.reshape(-1,1)
        model = LinearRegression()
        X = np.concatenate((X1, X2), axis=1)
        model.fit(X,VdCd)
        
        # Getting the estimated coefficients K and S
        K[j] = model.coef_[0]  # Coefficient for X1 conversion to mL/min
        S[j] = model.coef_[1]  # Coefficient for X2
        
        # Predicting VdCd using the fitted model
        predicted_VdCd = model.predict(X)
        
        # calculate predicted dialysate cocnetration
        for i in range(1,len(timesteps)):
            predicted_cd.loc[timesteps[i], solute] = (predicted_VdCd[i-1]+ V[timesteps[i-1]]*predicted_cd.loc[timesteps[i-1], solute])/V[timesteps[i]]
            
    new_index = pd.Index(np.arange(0, 241, 1))        
    predicted_cd = predicted_cd.reindex(new_index)
    predicted_cd = predicted_cd.interpolate(method = 'index', limit_direction = "both")
    predicted_cd.to_excel(f'predicted_cd/WM_F0p5_{ind}.xlsx')   
    #error calculation
    df2 = ((df_cd/cp-predicted_cd.loc[df_cd.index]/cp)**2)
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
    
    
    #calcualte RMSE
    obj_fn = sum(np.sqrt(df2.sum()/len(df_cd)))
    # "Use for plot"
    # fig, ax = plt.subplots(3,2, figsize = (12,18))
    # t = 240
    # #urea
    # df_cd['Urea'].plot( ax = ax[0,0], label = 'data', style = '.')
    # ax[0,0].plot(new_index,predicted_cd['Urea'], label = 'predicted')
    # # ax[0,0].text(0.6, 0.1, f'MTAC = {result["x"][0]:.2f} ml/min', transform=ax[0,0].transAxes)
    # ax[0,0].set_title("Urea")
    
    # #creatinine
    # df_cd['Creatinine'].plot( ax = ax[0,1], style = '.')
    # ax[0,1].plot(new_index,predicted_cd['Creatinine'])
    # # ax[0,1].text(0.6, 0.1, f'MTAC = {result["x"][1]:.2f} ml/min', transform=ax[0,1].transAxes)
    # ax[0,1].set_title("Creatinine")
    
    # #Sodium
    # df_cd['Sodium'].plot( ax = ax[1,0],  style = '.')
    # ax[1,0].plot(new_index,predicted_cd['Sodium'] )
    # # ax[1,0].text(0.6, 0.5, f'MTAC = {result["x"][2]:.2f} ml/min', transform=ax[1,0].transAxes)
    # ax[1,0].set_title("Sodium")
    
    # #Glucose
    # df_cd['Glucose'].plot( ax = ax[2,0], style = '.')
    # ax[2,0].plot(new_index,predicted_cd['Glucose'])
    # # ax[2,0].text(0.6, 0.5, f'MTAC = {result["x"][4]:.4f} ml/min', transform=ax[2,0].transAxes)
    # ax[2,0].set_title("Glucose")

    
    # fig.supxlabel("time, min")
    # fig.supylabel("Dialysate concentration, mmol")
    # plt.suptitle("Model 1 predictions of dialysate concentration")
    # plt.subplots_adjust(top=0.88,
    #                     bottom=0.11,
    #                     left=0.09,
    #                     right=0.9,
    #                     hspace=0.295,
    #                     wspace=0.215)
    et = time.time()
    comp_time = et-st
    row_data = [ind] + list(K) + list(S) + [obj_fn] + [comp_time]
    data.append(row_data)
   
cols = ['patient','MTAC_urea', 'MTAC_crea','MTAC_glu', 'MTAC_sodium',
        'sico_urea', 'sico_crea', 'sico_glu', 'sico_sodium', 'obj_fn', 'comp-time'] 
result = pd.DataFrame(data, columns=cols)
result.set_index('patient', drop = True, inplace = True)


with pd.ExcelWriter('fittedPS/fittedPS_WM_F0p5.xlsx', engine='openpyxl') as writer:
    result.to_excel(writer)