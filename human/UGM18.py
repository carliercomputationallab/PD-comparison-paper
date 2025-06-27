# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 16:23:23 2023
UNIFIED GRAFF MODEL all parameters
See unified_graff.py for explanation
@author: P70073624
"""


import multiprocessing
import numpy as np
import scipy
import pandas as pd
import time
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import Bounds


#%%
def calculate_penalty(x):
    #penalty for having sieving coefficient and fct beyond 0 and 1
    penalty = 0
    for param in x[4:12]:
        if param < 0:
            penalty += param ** 2
            
        elif param > 1:
            penalty += (param-1)**2
            
    return penalty
            
        

def objective(x, predicted_cd, Cp, V, df_cd):

    t = 240
    predicted_cd = rk(t, x, predicted_cd, Cp, V, df_cd)
    #error calculation
    df2 = ((df_cd/Cp-predicted_cd.loc[df_cd.index]/Cp)**2)
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
    penalty = calculate_penalty(x)
    return (sum(np.sqrt(df2.sum()/len(df_cd))) + penalty, predicted_cd)

#%%
def objective_fn(x, predicted_cd, Cp, V, df_cd):
    '''The objective function RMSE needed to be minimised'''
    # print(x)

    t = 240
    predicted_cd = rk(t, x, predicted_cd, Cp, V, df_cd)
    #error calculation
    df2 = ((df_cd/Cp-predicted_cd.loc[df_cd.index]/Cp)**2)
    df2['Glucose'] = ((df_cd['Glucose']/df_cd.loc[0, 'Glucose']-predicted_cd.loc[df_cd.index, 'Glucose']/df_cd.loc[0, 'Glucose'])**2)
    df2 = df2.astype('float64')
    penalty = calculate_penalty(x)
    print((df2.sum()/len(df_cd)))

    try:
        return sum(np.sqrt(df2.sum()/len(df_cd))) + penalty #RMSE
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


# the differential equations for Graff model
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
    L = 0.3 # lymphatic flow rate ml/min
    QU = V[t+1] - V[t] + L #UF rate ml/min delta t is 1
    beta = QU * SiCo/MTAC.ravel() #ml/min
    f = np.array([0 if b > 30 else 0.5 if b == 0 else (1/b)-1/(np.exp(b)-1) for b in beta]) #check model 2 for explanation
    MC = cp-f*(cp-cd)
    Cl = cp if L < 0 else cd
    dxdt = (MTAC * (fct * cp - cd) + SiCo * QU * MC - L * Cl - cd * QU)/V[t].ravel()
    
    
    return dxdt

def multiprocessing_func(pfile):
    
    st = time.time()
    Nx = 12  #4MTACs, 4Sicos, 4 fcts
    Cp = data[pfile]['cp']
    V = data[pfile]['V']
    df_cd = data[pfile]['df_cd']
    predicted_cd = pd.DataFrame(columns= ["Urea", "Creatinine", "Glucose", "Sodium"], dtype = float)
    predicted_cd.loc[0]=df_cd.loc[0]
    optimised_values = np.empty(Nx)
    obj_fn = []
    print(pfile)
    
    for var in range(10):
        
        x0 = np.random.randint(1, 21, size=Nx)
        lbound = [0.01, 0.01, 0.01, 0.01,  0, 0,  0, 0, 0, 0,  0, 0]
        ubound = [200, 200, 200, 200, 1, 1, 1, 1, 1, 1, 1, 1] 
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
#%%    
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
  


#%%
'''
single patient simulation
'''
# pfile = data['PDPC_754']
# Nx = 3 #for MTACs
# Cp = pfile['cp']
# V = pfile['V']
# df_cd = pfile['df_cd']
# Vr = pfile['Vr']
# V_fill = pfile['Vfill']
# LS = 0.02


# predicted_cd = pd.DataFrame(columns= ["Urea", "Creatinine", "Glucose"], dtype = float)

# predicted_cd.loc[0]=df_cd.loc[0]

# predicted_V = np.zeros(241)
# predicted_V[0] = V[0]

# optimised_values = np.empty(Nx)
# obj_fn = []

# for var in range(10):
#     print(var)
#     #Define initial initial_guess
#     x0 = np.array(random.sample(range(1, 50), Nx))#np.concatenate([np.array(random.sample(range(1, 50), 3)), np.array(random.sample(range(1, 46), 1))/50])
    
#     '''SLSQP optimisation'''
#     result = scipy.optimize.minimize(objective_fn, x0, args = (predicted_cd, Cp, V, df_cd, Vr, V_fill, predicted_V, LS),
#             method='SLSQP', bounds = [(0.01, 200), (0.01, 200),(0.01, 200) ],
#             options = {"maxiter" : 1000, "disp": True})
    
#     #gather all optimised values
#     optimised_values = np.vstack((optimised_values,result['x'].tolist()))
    
    
    
#     #get total error
#     obj_fn.append(result['fun'])
#%%

if __name__ == '__main__':

    pool = multiprocessing.Pool(9)
    result_list = pool.map(multiprocessing_func,patientlist)
    pool.close()

    OF = []
    OV = np.empty(len(patientlist),dtype=object)
    UF = np.empty(len(patientlist),dtype=object)
    comp_time = np.empty(len(patientlist),dtype=float)
    Cerr = np.empty(len(patientlist),dtype=float)
    
    for i in range(len(result_list)):
        OF.append(min(result_list[i][1]))
        OV[i] = result_list[i][0][np.argmin(result_list[i][1])+1]
        comp_time[i] = result_list[i][2]
        pfile = patientlist[i]
        
        Cp = data[pfile]['cp']
        V = data[pfile]['V']
        df_cd = data[pfile]['df_cd']
        predicted_cd = pd.DataFrame(columns= ["Urea", "Creatinine", "Glucose", "Sodium"], dtype = float)
        predicted_cd.loc[0]=df_cd.loc[0]


        
        _, predicted_cd = objective(OV[i], predicted_cd, Cp, V, df_cd) #collect all volume flow through pores
        predicted_cd.to_excel(f'predicted_cd/UGM18_{pfile}.xlsx') 
        # fig, ax = plt.subplots(1,4, figsize = (18,6))
        # t = 240
        # #urea
        # df_cd['Urea'].plot( ax = ax[0], label = 'data', style = '.')
        # ax[0].plot(np.arange(t+1),predicted_cd['Urea'], label = 'predicted')
        # # ax[0,0].text(0.6, 0.1, f'MTAC = {result["x"][0]:.2f} ml/min', transform=ax[0,0].transAxes)
        # ax[0].set_title("Urea")
        
        # #creatinine
        # df_cd['Creatinine'].plot( ax = ax[1], style = '.')
        # ax[1].plot(np.arange(t+1),predicted_cd['Creatinine'])
        # # ax[0,1].text(0.6, 0.1, f'MTAC = {result["x"][1]:.2f} ml/min', transform=ax[0,1].transAxes)
        # ax[1].set_title("Creatinine")
        
       
        
        # #Glucose
        # df_cd['Glucose'].plot( ax = ax[2], style = '.')
        # ax[2].plot(np.arange(t+1),predicted_cd['Glucose'])
        # # ax[2,0].text(0.6, 0.5, f'MTAC = {result["x"][4]:.4f} ml/min', transform=ax[2,0].transAxes)
        # ax[2].set_title("Glucose")
        
        # ax[3].plot(np.arange(t+1), pfile['V'], lw = 3, label = 'Expt')
        # ax[3].plot(np.arange(t+1), UF[i], label = 'Sim')
        # ax[3].set_ylabel('Volume, ml')
        # ax[3].legend()
        
        # fig.supxlabel("time, min")
        # fig.supylabel("Dialysate concentration, mmol")
        # plt.suptitle(datafile.index[i])
        # plt.savefig(datafile.index[i]+'.png')
        
        # with pd.ExcelWriter('fittedPS/predicted_cd_UGM18_PDPC_754.xlsx', engine='openpyxl') as writer:
        #     predicted_cd.to_excel(writer) 
        
        
        
    #%% 
    # Data analysis and results processing
    # ...
    # The code collects and processes simulation results
    # It calculates various metrics and saves them to CSV files
    cols = ['MTAC_urea', 'MTAC_crea','MTAC_glu', 'MTAC_sodium',
            'fct_urea', 'fct_crea','fct_glu', 'fct_sodium',
            'sico_urea', 'sico_crea', 'sico_glu', 'sico_sodium'] 
    PS = pd.DataFrame([arr for arr in OV], columns=cols)  
    PS['patient'] = patientlist
    PS['obj_fn'] = OF
    PS['comp-time'] = comp_time
    PS.set_index('patient', drop = True, inplace = True)
    #PS.to_csv('fittedPS/fittedPS_unifiedGraff_human_18p_PDPC_754_1.csv')
    with pd.ExcelWriter('fittedPS/fittedPS_UGM18.xlsx', engine='openpyxl') as writer:
        PS.to_excel(writer) 