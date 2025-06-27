# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 15:59:57 2022
New Waniewski model fitting
Check corresponding file in human for details
the only difference is 6 solutes instead of 3.
@author: P70073624
"""

import numpy as np
import os
import pandas as pd
from fnmatch import fnmatch
from scipy.interpolate import interp1d
from sklearn.preprocessing import MinMaxScaler
import time
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.integrate import quad

def integrand1(t, df):
    t = int(t)   
    return Cp[t,j] - df.loc[t, solute]

def integrand2(t, df):
    t = int(t)
    return Qv* ((1- F)* Cp[t,j]+ F* df.loc[t, solute])
    

"Get all MTACs for the pig in question"
root = 'patient_files/'
pattern = "*.csv"
patientlist = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            patientlist.append(os.path.join(path, name))
  
# patientlist = ["patient_files/pig4\session2\P29.csv"]
scaler = MinMaxScaler()
obj_fn = []

"Get MTAC"       
solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]

patient_no = []
MTAC_G = pd.Series(dtype=float) 
time_Garred = np.empty(len(patientlist))
data = []
for k,pfile in enumerate(patientlist):  
    st = time.time()      

    p = pfile.split("\\")[-1] #to get the file name
    print(p)
    patient_no.append(p)
    df = pd.read_csv(pfile,skiprows = range(0,16), delimiter = "," \
                     , encoding= 'unicode_escape')
    # print(df.head())
    '''Plasma solute concentration'''
    df_cp = pd.DataFrame(index=[0, 120, 240],columns= solutes, dtype = float)
    c_prot = np.nanmean(df['Total protein'].iloc[1:4].astype('float64').to_numpy()) #in g/L
    
    df_cp = df[solutes].iloc[1:4].copy()#blood plasma concentration

    if pd.isna(df_cp["Potassium"]).all():
        if float(df["Potassium"].iloc[4]):
            df_cp.loc[:,"Potassium"]=float(df["Potassium"].iloc[4])*0.96
        else:
            df_cp.loc[:,"Potassium"]=0.96*4.0 # 0.96 is Donnan correction       
        
    for column in df_cp:
        df_cp[column] = df_cp[column].astype('float64')*1/(0.984 - 0.000718*c_prot) #plasma water correction, c_prot should be in g/L
    index = pd.Series([0,120,240])
    df_cp = df_cp.set_index([index])
    df_cp = df_cp.interpolate(method = 'index', limit_direction = "both")  
    df_cp.loc[:,"Sodium"] *=0.96
    df_cp.loc[:, "Creatinine"]  *= 0.001
    # uses the interpolation using the values of the indices and interpolates in both directions
    f_cp = interp1d(df_cp.index, df_cp, axis = 0)
    interpolated_cp = f_cp(range(0,241))
    Cp = interpolated_cp
    

    '''dialysate solute concentration'''
    df_cd = pd.DataFrame(columns = solutes,dtype = float)
    df_cd = df[solutes].iloc[10:18].copy()  #dialysate concentration
    for column in df_cd:
        df_cd[column] = df_cd[column].astype('float64')  
    df_cd.loc[:, "Creatinine"]  *= 0.001   
    index = pd.Series([0,10,20,30,60,120,180, 240])
    df_cd = df_cd.set_index([index])
    df_cd = df_cd.interpolate(method = 'index', limit_direction = "both")
    # Create a new index from 0 to 240 with a step of 1
    new_index = pd.Index(np.arange(0, 241, 1))

    df_cd2 = df_cd.copy(deep = True)
    # Reindex the DataFrame
    df_cd2 = df_cd2.reindex(new_index)
    df_cd2 = df_cd2.interpolate(method = 'index', limit_direction = "both")
    
    '''dialysate volume'''
    V = pd.read_csv(pfile,skiprows = range(0,45), delimiter = ",", \
                        encoding= 'unicode_escape')[["IP volume T=0 (mL)","IP volume T=240 (mL)"]].iloc[0] # IPV measured from haemoglobin
    # print(df_V)
    df_V = pd.read_csv(pfile,skiprows = range(0,60), delimiter = ",", \
                        encoding= 'unicode_escape')[["Time","IPV"]].iloc[2:10]
        
    
    df_V['IPV'] = df_V['IPV'].replace(["#REF!",'#DIV/0!', '#VALUE!'], np.nan).interpolate()
    df_V =df_V.astype('float64')
    df_V.set_index('Time', inplace = True)
    if np.isnan(df_V.iloc[0].iloc[0]):
        df_V.loc[0] = float(V.iloc[0])
        
    if np.isnan(df_V.iloc[-1].iloc[0]):
        df_V.iloc[-1] = float(V.iloc[1])
        
    df_V['IPV'] = df_V['IPV'].interpolate()


    #Linear interpolation to find values of V at all times
    f_V = interp1d(df_V.index, df_V, axis = 0)
    interpolated_V = f_V(range(0,241))
    V = interpolated_V.flatten()
    
    # print(df_V)
    
    F = 1
    timesteps = np.array(df_cd.index)
    X1 = np.empty(len(timesteps)-1)
    X2 = np.empty(len(timesteps)-1)
    VdCd = np.empty(len(timesteps)-1)
    K = np.empty(len(solutes))
    S = np.empty(len(solutes))
    predicted_cd = pd.DataFrame(columns= ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"], dtype = float)
    
    predicted_cd.loc[0]=df_cd.loc[0] #the value of predicted cd0 comes from the initial guess
    for j, solute in enumerate(solutes):
        for i in range(1,len(timesteps)):
            X1[i-1] = quad(integrand1, timesteps[i-1], timesteps[i], df_cd2)[0]
            # X1[i-1] = quad(integrand1, timesteps[i-1], timesteps[i], limit = 100)[0]
            # Qv = (V[timesteps[i]] - V[timesteps[i-1]])/(timesteps[i]-timesteps[i-1]) + 0.0007
            Qv = (V[240]-V[0])/240 - 0.3 *240
            # print(Qv)
            # X2[i-1] = Qv* (timesteps[i]-timesteps[i-1]) * ((1-F)*Cp[timesteps[i],j]+F*df_cd.loc[timesteps[i], solute] - ((1-F)*Cp[timesteps[0],j]+F*df_cd.loc[timesteps[0], solute]))
            X2[i-1] = quad(integrand2, timesteps[i-1], timesteps[i], df_cd2)[0]
            VdCd[i-1] = V[timesteps[i]]*df_cd.loc[timesteps[i], solute] - V[timesteps[i-1]]*df_cd.loc[timesteps[i-1], solute]
        
        X1 = X1.reshape(-1,1)
        # plt.figure(1)
        # plt.plot(X1, VdCd)
        X2 = X2.reshape(-1,1)
        model = LinearRegression()
        X = np.concatenate((X1, X2), axis=1)
        model.fit(X,VdCd)
        # Getting the estimated coefficients K and S
        K[j] = model.coef_[0]  # Coefficient for X1 conversion to mL/min
        S[j] = model.coef_[1]  # Coefficient for X2

        # Plotting the fits for each patient
        # plt.figure(figsize=(10, 6))
        # plt.scatter(timesteps[1:], VdCd, label='Actual Data', color='blue')
        
        # Predicting VdCd using the fitted model
        predicted_VdCd = model.predict(X)
        
        # plt.plot(timesteps[1:], predicted_VdCd, label='Fitted Data', color='red')
        # plt.xlabel('Timesteps')
        # plt.ylabel('VdCd')
        # plt.title(f'Fit for Patient {p}, solute{solute}')
        # plt.legend()
        # plt.show()
        
        for i in range(1,len(timesteps)):
        
            predicted_cd.loc[timesteps[i], solute] = (predicted_VdCd[i-1]+ V[timesteps[i-1]]*predicted_cd.loc[timesteps[i-1], solute])/V[timesteps[i]]
            
    predicted_cd = predicted_cd.reindex(new_index)
    predicted_cd = predicted_cd.interpolate(method = 'index', limit_direction = "both")
    predicted_cd.to_excel(f'predicted_cd/WM_F1_{pfile[28:]}.xlsx')
    df2 = ((df_cd/df_cp.loc[0]-predicted_cd.loc[df_cd.index]/df_cp.loc[0])**2).astype('float64')
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')

    obj_fn = sum(np.sqrt(df2.sum()/len(df_cd)))
    # "Use for plot"
    # fig, ax = plt.subplots(3,2, figsize = (12,18))
    # t = 240
    # #urea
    # df_cd['Urea'].plot( ax = ax[0,0], label = 'data', style = '.')
    # ax[0,0].plot(timesteps,predicted_cd['Urea'], label = 'predicted')
    # # ax[0,0].text(0.6, 0.1, f'MTAC = {result["x"][0]:.2f} ml/min', transform=ax[0,0].transAxes)
    # ax[0,0].set_title("Urea")
    
    # #creatinine
    # df_cd['Creatinine'].plot( ax = ax[0,1], style = '.')
    # ax[0,1].plot(timesteps,predicted_cd['Creatinine'])
    # # ax[0,1].text(0.6, 0.1, f'MTAC = {result["x"][1]:.2f} ml/min', transform=ax[0,1].transAxes)
    # ax[0,1].set_title("Creatinine")
    
    # #Sodium
    # df_cd['Sodium'].plot( ax = ax[1,0],  style = '.')
    # ax[1,0].plot(timesteps,predicted_cd['Sodium'] )
    # # ax[1,0].text(0.6, 0.5, f'MTAC = {result["x"][2]:.2f} ml/min', transform=ax[1,0].transAxes)
    # ax[1,0].set_title("Sodium")
    
    # #Phosphate
    # df_cd['Phosphate'].plot( ax = ax[1,1], style = '.')
    # ax[1,1].plot(timesteps,predicted_cd['Phosphate'] )
    # # ax[1,1].text(0.6, 0.1, f'MTAC = {result["x"][3]:.2f} ml/min', transform=ax[1,1].transAxes)
    # ax[1,1].set_title("Phosphate")
    
    # #Glucose
    # df_cd['Glucose'].plot( ax = ax[2,0], style = '.')
    # ax[2,0].plot(timesteps,predicted_cd['Glucose'])
    # # ax[2,0].text(0.6, 0.5, f'MTAC = {result["x"][4]:.4f} ml/min', transform=ax[2,0].transAxes)
    # ax[2,0].set_title("Glucose")
    
    # #Potassium
    # df_cd['Potassium'].plot( ax = ax[2,1], style = '.')
    # ax[2,1].plot(timesteps,predicted_cd['Potassium'])
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
    et = time.time()
    time_Garred[k] = et-st
    row_data = [p] + list(K) + list(S) + [obj_fn] + [time_Garred[k]]
    data.append(row_data)
    
    

    
    

    
    
cols = ['patient','MTAC_urea', 'MTAC_crea','MTAC_sodium', 'MTAC_phosphate','MTAC_glu', 'MTAC_potassium',
        'sico_urea', 'sico_crea','sico_sodium', 'sico_phosphate', 'sico_glu', 'sico_potassium', 'obj_fn', 'comp time'] 
result = pd.DataFrame(data, columns=cols)
result.set_index('patient', drop = True, inplace = True)
"Transpose the dataframe and remove the first empty line and index it with patient name"

with pd.ExcelWriter('fittedPS/fittedPS_WM_F1.xlsx', engine='openpyxl') as writer:
    result.to_excel(writer)
