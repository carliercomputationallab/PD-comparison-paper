# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:27:18 2022
Extract cp, cd and Vfill and Vres from patient files
@author: P70073624
"""
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np

def input_values(pfile):
    
    solutes = ["Urea", "Creatinine", "Sodium", "Phosphate", "Glucose", "Potassium"]
    t = 240 #min
    #p = input("Enter patient's file number")
    #pfile = p +".csv"
    #pfile = "2019415_1.csv"
    p = pfile.split(".")[0]
    # print(p)
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
    # print(df_cp)

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
    #Linear interpolation to find values of cp at all times
    
    f_cp = interp1d(df_cp.index, df_cp, axis = 0)
    interpolated_cp = f_cp(range(0,t+1))

    #Linear interpolation to find values of V at all times
    f_V = interp1d(df_V.index, df_V, axis = 0)
    interpolated_V = f_V(range(0,t+1))

    #predicted_cd[0] = df_cd["cd"][0]
    
    Cp = interpolated_cp
    V = interpolated_V.flatten()
    
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
    
    return predicted_cd, Cp, V, df_cd, Vr, V_fill, df_V

