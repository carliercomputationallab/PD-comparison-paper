****directory****
patient files -> has all the data separated by pig and session
spider_data -> all fitted parameters for teh three iterations of each model 1- 9 and both data splits 7:4 and 6:5

****Python files****
TPM = three-pore model
GM = Garred model
SWM = simplified Waniewski model
WM = Waniewski Model
UGM = unified Graff model (best fit for fct and SiCo)
UGM18 = unified GRaff model (all random parameters)
 -> has parallel processing routine for all model
-> save PS, UF and computational time for the particular model

sse_* = calculate solute-specific error for all models


values - > this is where we extract all plasma and dialysate concentration directly from the patient files

UF values -> *for model 7 V fitting only* In addition to all things in values, we also extract delP, Vres, L

compare models AIC -> different ways of comparing the models
by RMSE error, AIC and fitted parameter values.

compare SSE -> compare solute specific error

