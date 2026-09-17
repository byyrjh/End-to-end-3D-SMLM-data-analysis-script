clear;clc;
load('mu_vec_SM0.mat')
load('seg_config_RMSE_optimal_SM0.mat')
load('fitting_result_SM_m0_round2.mat')
% config 1 z = -800:800
% config 2 z = -400 0 400 800
% config 3 z = -800 -400 0 400 
% When number of segment slices = 4 origin in z is fixed at slice 3
% thus in config 2 mu = mu - 400;
mu_vec(seg_config_RMSE == 2) = mu_vec(seg_config_RMSE == 2) -400;
fit_err = fitting_results(3,:)*10 + fitting_results(6,:)*10 - mu_vec';
figure
histogram(fit_err)