Figure 2: compare RMSE
Run compare_RMSE_pigs_human.py

Figure 3: compare solute specific error
Run compare_SSE_pigs_human.py
Prior to this run,
sse_m7, sse_m8, sse_m9, sse_unifiedGraff, sse_unifiedGraff_18p, sse_Waniewski

Figure 4: Specific pig concentration profiles


Figure 5: Specific human concentration profiles
human/specific_human/plot_predicted_cd_specific.py

Figure 6: Compare fitted parameters
Run compare_fitted_parameters.py

Figure 7: Treament Score
human/specific_human/optimise_patient_clearance.py

Figure 8: Correlation matrix
human/correlation_matrix.py

Folder: human
model7-train-parallel.py -- three-pore model (TPM)
model8.py -- Garred model (GM)
model9.py -- simple Waniewski model (SWM)
new_waniewski_model.py -- Waniewski model (WM)
unified_graff -- Graff model with the best model parameters (SiCo and fct) of each solute  (UGM)
unified_graff_18p -- Graff model with all unknown parameters (UGM-18)
compare_models_AIC -- comparison of fitted parameters and computational time; write mean and SD of the fitted parameters of all 20 human
correlation_matrix.py -- human demographics and correlation with fitted parameters
sse_xx -- calculate solue specific error with the predicted MTACS (where xx is m7 (TPM), m8 (Garred), m9 (SWM), unified Graff, UGM18, WM)

Folder: human/specific_human
xx_specific - fit dialysate concetraiton profiles for only 1 patient
compare_fittedPS_per_patient -- compile fitted parameters from 6 models
model7-Vfitting-withLS -- Fitting LS to predict volume profile (procedure in Figure S1)
optimise_patient_clearance -- Figure 7 (heat map treatment score)
plot_predicted_cd -- predict human dialysate concentration profiles Figure supplementary S7
plot_predicted_cd_specific -- predict human dialysate concentration profiles Figure 5