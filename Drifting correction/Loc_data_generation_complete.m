clear;clc;close all;
setting_para.smooth_f = 20;% determined by low vibrational frequency  % 30  % 6
setting_para.exposuretime = 20;  %ms  % 30  % 20   
setting_para.backlash_thres = 1500; %nm  distance between FM in m0 and m1 image stacks  % 1000  % 1500
setting_para.hyper_piece_idx = 1;
setting_para.xybinsize = 1000; % unit nm
correct_SM = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  get correction coef  hyperstack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  segment size  
FM_outlier = cell(2,1);  % inspect fiducial marker traces in FM_loc_prec.m and kickout outliers
FM_outlier{1,1} = [2 5 13 14 15 17 18];  
FM_outlier{2,1} = [6 14];
setting_para.piecepos = [1 13 16 18 19 58 59 60]; % including start == 1 and end == (number of hyperstack +1) each index have jump at its initial volume
data_path='N:\Lucas\10_LLSM_CudaFitter\03_HeLa_SMLM\2023-08-16-Hela-aTUB\2023-08-16_15-41-01_Hela-ATUB-230815-AB-FOV7-002\';  
setting_para.stationary_pos=0;  % unstable volumes are kicked out before fitting stage
setting_para.outlier_thres_xy = 120; % kick out spikes in FM trajectory
setting_para.outlier_thres_z = 120; % kick out spikes in FM trajectory
setting_para.pixel_size_z = 10;  % nm
setting_para.pixel_size_xy = 100;  % nm
setting_para.num_photon_factor=30;
setting_para.FM_outlier = FM_outlier; % kick out FM outliers
setting_para.origin_drifting=1; % define origin of FM trajectory, which respect to which the SM data is corrected
setting_para.data_path = data_path;
setting_para.coef_vis = true;
setting_para.diff_vis = true;
setting_para.close_loop = false;  % if sample piezo stage is controlled by close loop (if true, traces are smoothed in segment)
setting_para.m1 = true;
setting_para.SlicesPerVol = 50;
%%
load(strcat(data_path,'PSFdata_FM_m0.mat')); % M:\Hao\2022\2022_11_18\cell1\
setting_para.imgsize = size(data.data{1,1}.xyMIP,1);
setting_para.slice_per_stack = data.slicepS;
setting_para.step_size = data.stepsize*1000;  % nm
setting_para.vol_per_hyper = length(data.data)/length(data.path);
setting_para.num_hyper = length(data.path);
datasave_path = strcat(data_path,'segment_data\');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     diff_x and diff_y are piecewise funtions      %%%%%%%%%%%%%
%%%     diff_x_1 = a*t + b
%%%     diff_x_1 = (coef_ax_1(1)*z + coef_ax_1(2))*t+(coef_bx_1(1)*z^2+coef_bx_1(2)*z+coef_bx_1(3))
%%%     diff_z = coef_z(1)*t + coef_z(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     update version in 5.3.2024
%%%     diff_x and diff_y are corrected for each time t by quadratic curve
%%%     fitting
%%%     diff_z are corrected using mean z trajectories after light sheet plane fluctuation correction
%%%     delta ls_os = apparent ls_os - light sheet plane fitting
% FM_trace_m0 = get_FM_trace(datasave_path, 'm0', setting_para);
SM_data_m1_cor.x = [];
SM_data_m1_cor.y = [];
SM_data_m1_cor.z = [];
SM_data_m1_cor.inten = [];
SM_data_m1_cor.t = [];

if setting_para.m1
    [coef FM_trace_m0]= get_correction_coef(datasave_path,setting_para);    
end
if correct_SM
    SM_data_m0 = get_SM_data(datasave_path, 'm0', setting_para,coef);
    SM_data_m0_cor = drifting_correct(SM_data_m0,FM_trace_m0,setting_para, 'm0');
    SM_data_m1 = get_SM_data(datasave_path, 'm1', setting_para,coef);
    SM_data_m1_cor = drifting_correct(SM_data_m1,FM_trace_m0,setting_para, 'm1');
end
%%
theta = 33/180*pi;
t = [ 1          0            0;
    0      cos(theta)  -sin(theta);
    0      sin(theta)   cos(theta) ];
x_pos_full = [SM_data_m0_cor.x; SM_data_m1_cor.x];
y_pos_full = [SM_data_m0_cor.y; SM_data_m1_cor.y];
z_pos_full = [SM_data_m0_cor.z; SM_data_m1_cor.z];
inten_full = [SM_data_m0_cor.inten; SM_data_m1_cor.inten];
vol_idx_full = [SM_data_m0_cor.t; SM_data_m1_cor.t];
full_SM_data=[x_pos_full y_pos_full z_pos_full inten_full vol_idx_full];
full_SM_data(:,4) = 1;
full_SM_data_bio = full_SM_data;
full_SM_data_bio(:,1:3)=(t*full_SM_data(:,1:3)')';

tg=fopen(strcat(datasave_path,'LS_mechanical_cor_bidir.txt'),'wt');
for j=1:size(full_SM_data_bio,1)
    fprintf(tg,'%f\t%f\t%f\t%f\t%f\n',full_SM_data_bio(j,1),full_SM_data_bio(j,2),full_SM_data_bio(j,3),full_SM_data_bio(j,4),full_SM_data_bio(j,5));
end
fclose(tg);
