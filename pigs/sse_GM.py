# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 16:08:24 2023

@author: P70073624
"""



import numpy as np
import os
import pandas as pd
from fnmatch import fnmatch
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

"Get all MTACs for the pig in question"
root = 'patient_files/'
pattern = "*.csv"
patientlist = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            patientlist.append(os.path.join(path, name))

# patientlist = ['patient_files/pig2/session1/P10.csv']
solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
mtac = pd.read_excel('fittedPS/fittedPS_GM.xlsx', index_col = 0)
sse = pd.DataFrame(columns=solutes)
"Get MTAC"       

# patientlist= ['P9.csv', 'P15.csv', 'P11.csv', 'P21.csv', 'P10.csv', 'P23.csv', 'P8.csv']

for i,pfile in enumerate(patientlist):  
  
    t = 240 #min
    p = pfile[28:] #to get the file name
    
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
    print(df_cp.mean())
    

    '''dialysate solute concentration'''
    df_cd = pd.DataFrame(columns = solutes,dtype = float)
    df_cd = df[solutes].iloc[10:18].copy()  #dialysate concentration
    for column in df_cd:
        df_cd[column] = df_cd[column].astype('float64')  
    df_cd.loc[:, "Creatinine"]  *= 0.001   
    index = pd.Series([0,10,20,30,60,120,180, 240])
    df_cd = df_cd.set_index([index])
    #using .values here to copy only the values, otherwise it tries to match the indices of df and df_cd and it doesnt work
    df_cd = df_cd.interpolate(method = 'index', limit_direction = "both")
    # print(df_cd)

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
    interpolated_V = f_V(range(0,t+1))
    V = interpolated_V.flatten()
    
    x = mtac.loc[p][0:6]
    
    x.index = df_cp.columns
    
    predicted_cd = pd.DataFrame(columns= solutes, dtype = float)
    predicted_cd.loc[0]=df_cd.loc[0]
    
    
    for t in range(1,241):
        predicted_cd.loc[t] = df_cp.mean()-(df_cp.mean()-predicted_cd.loc[0])*(V[0]/V[t])*np.exp(-x/V.mean()*t)

        
    df2 = ((df_cd/df_cp.loc[0]-predicted_cd.loc[df_cd.index]/df_cp.loc[0])**2).astype('float64')
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
  
    err = np.sqrt(df2.sum()/len(df_cd))

    sse.loc[p] = err
    
with pd.ExcelWriter('fittedPS/sse_GM.xlsx', engine='openpyxl') as writer:
    sse.to_excel(writer)   

