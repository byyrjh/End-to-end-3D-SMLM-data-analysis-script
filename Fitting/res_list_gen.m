%% warning!!!!!
% this code should be run for only once!!!
% map_ptr_z_SM_m0.mat already updated
clear;clc;
load(strcat('map_ptr_z_SM_m0.mat'));
load(strcat('N:\Lucas\10_LLSM_CudaFitter\03_HeLa_SMLM\2023-08-16-Hela-aTUB\2023-08-16_15-41-01_Hela-ATUB-230815-AB-FOV7-002\segment_data\redefined RMSE\reduced roi 11_11\seg_config_RMSE_optimal_SM0.mat'));
idx_case2 = find(seg_config_RMSE==2);
map_ptr_z_SM(idx_case2) = map_ptr_z_SM(idx_case2)+1;
save('map_ptr_z_SM_m0.mat','map_ptr_z_SM');

load(strcat('map_ptr_z_SM_m1.mat'));
load(strcat('N:\Lucas\10_LLSM_CudaFitter\03_HeLa_SMLM\2023-08-16-Hela-aTUB\2023-08-16_15-41-01_Hela-ATUB-230815-AB-FOV7-002\segment_data\redefined RMSE\reduced roi 11_11\seg_config_RMSE_optimal_SM1.mat'));
idx_case2 = find(seg_config_RMSE==2);
map_ptr_z_SM(idx_case2) = map_ptr_z_SM(idx_case2)+1;
save('map_ptr_z_SM_m1.mat','map_ptr_z_SM');