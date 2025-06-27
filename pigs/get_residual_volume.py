# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 15:26:04 2024

@author: P70073624
"""

import pandas as pd
import numpy as np
import os
from fnmatch import fnmatch
from scipy.interpolate import interp1d

root = 'patient_files/'
pattern = "*.csv"
patientlist = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            patientlist.append(os.path.join(path, name))

for pfile in patientlist:

    solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
    t = 240 #min
    #p = input("Enter patient's file number")
    #pfile = p +".csv"
    #pfile = "2019415_1.csv"
    p = pfile.split(".")[0]
    # print(p)
    df = pd.read_csv(pfile,skiprows = range(0,15), delimiter = "," \
                     , encoding= 'unicode_escape')
    
    
    '''dialysate solute concentration'''
    df_cd = pd.DataFrame(columns = solutes,dtype = float)
    df_cd = df[solutes].iloc[8:19].copy()  #dialysate concentration
    for column in df_cd:
        df_cd[column] = df_cd[column].astype('float64')  
    df_cd.loc[:, "Creatinine"]  *= 0.001   
    index = pd.Series([0,10,20,30,60,120,180, 240])
    df_cd = df_cd.set_index([index])
    #using .values here to copy only the values, otherwise it tries to match the indices of df and df_cd and it doesnt work
    df_cd = df_cd.interpolate(method = 'index', limit_direction = "both")
    # print(df_cd)
    
    df_cd['Dextran'] = df.iloc[0, 8:19]
    
    '''dialysate volume'''
    df_V = pd.read_csv(pfile,skiprows = range(0,45), delimiter = ",", \
                       encoding= 'unicode_escape')[["IP volume T=0 (mL)","IP volume T=240 (mL)"]].iloc[0] # IPV measured from haemoglobin

    
    #Linear interpolation to find values of V at all times
    f_V = interp1d([0,240], df_V)
    interpolated_V = f_V(range(0,t+1))
    
    #predicted_cd[0] = df_cd["cd"][0]
    
    V = interpolated_V/1000
    
    cd0_m, MTAC_m, fct_m, SiCo_m, QL_m, AIC =[np.empty(6, dtype = object) for _ in range(6)]
    
    #cols = ["model", "cd0", "MTAC","fct","SiCo","QL", "Guess no.", "evaluation", "AIC_score", "Initital guess"]
    
    predicted_cd = pd.DataFrame(columns= ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"], dtype = float)
    
    predicted_cd.loc[0]=df_cd.loc[0] #the value of predicted cd0 comes from the initial guess
    
    # Residual volume ... total protein
    Vr = float(pd.read_csv(pfile,skiprows = range(0,45), delimiter = ",", \
                       encoding= 'unicode_escape')["RV before SPA (mL)"].iloc[1]) #ml # Changing back to normal values
    # fill volume
    V_fill = float(pd.read_csv(pfile,skiprows = range(0,39), delimiter = ",", \
                       encoding= 'unicode_escape')["Volume instilled (mL)"].iloc[1]) #ml Changing back to normal values
    if np.isnan(V_fill):
        V_fill = 2000 # mL. in case there are no measurement values
    df_drain = pd.DataFrame(columns = ["Creatinine", "Total protein", "Albumin"],dtype = float)
    df_drain = df[["Creatinine", "Total protein", "Albumin"]].iloc[8:11].copy()  #dialysate concentration
    for column in df_drain:
        df_drain[column] = df_drain[column].astype('float64') 
    df_drain.loc[:, "Creatinine"]  *= 0.001
    
    for col in df_drain.columns:   
        if np.isnan(Vr):
            Vr = V_fill*df_drain.loc[10,col]/(df_drain.loc[10,col]-df_drain.loc[8,col])
    
    if np.isnan(Vr):
        Vr = 200 # mL. in case there are no measurement values