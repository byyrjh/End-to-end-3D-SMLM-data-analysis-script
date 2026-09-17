clear;clc;
fitting_round = 1;
fit_offset = 1;
scan_mode = 'm0';
scriptFolder = fileparts(mfilename('fullpath'));
packageFolder = fileparts(scriptFolder);
data_path = fullfile(packageFolder,'data\segment_data');
load(strcat(data_path,'\seg_config_RMSE_optimal_SM_',scan_mode,'.mat'));
load(strcat(data_path,'\seg_data_SM_',scan_mode,'.mat'));
load(strcat(data_path,'\map_ptr_x_SM_',scan_mode,'.mat'));
load(strcat(data_path,'\map_ptr_y_SM_',scan_mode,'.mat'));
if fitting_round>1
    load(strcat(data_path,'\fitting_init_SM_',scan_mode,'.mat'));
end
seg_5_idx = find(seg_config_RMSE == 1);
seg_4down_idx = find(seg_config_RMSE == 2);
seg_4up_idx = find(seg_config_RMSE == 3);
ptr_x_5 = map_ptr_x_SM(seg_5_idx);
ptr_y_5 = map_ptr_y_SM(seg_5_idx);
ptr_x_4_down = map_ptr_x_SM(seg_4down_idx);
ptr_y_4_down = map_ptr_y_SM(seg_4down_idx);
ptr_x_4_up = map_ptr_x_SM(seg_4up_idx);
ptr_y_4_up = map_ptr_y_SM(seg_4up_idx);
[size_xy, ~, ~] = size(seg_data_SM);
seg_5 = zeros(size_xy,size_xy,5*length(seg_5_idx));
seg_4_down = zeros(size_xy,size_xy,4*length(seg_4down_idx));
seg_4_up = zeros(size_xy,size_xy,4*length(seg_4up_idx));
for i = 1:length(seg_5_idx)
    seg_temp = seg_data_SM(:,:,(seg_5_idx(i)-1)*5+1:seg_5_idx(i)*5);
    seg_5(:,:,(i-1)*5+1:i*5) = seg_temp;
end
for i = 1:length(seg_4down_idx)
    seg_temp = seg_data_SM(:,:,(seg_4down_idx(i)-1)*5+1:seg_4down_idx(i)*5);
    seg_4_down(:,:,(i-1)*4+1:i*4) = seg_temp(:,:,2:5);
end
for i = 1:length(seg_4up_idx)
    seg_temp = seg_data_SM(:,:,(seg_4up_idx(i)-1)*5+1:seg_4up_idx(i)*5);
    seg_4_up(:,:,(i-1)*4+1:i*4) = seg_temp(:,:,1:4);
end
%% fitting parameter setting
fitting_results = zeros(6,length(map_ptr_x_SM));
data_dir = data_path;
cali_dir = fullfile(packageFolder,'data\');
load(strcat(cali_dir,'setup_calibration\camera1_cali.mat'));
load(strcat(cali_dir,'coeff_det_exper1.mat'));
load(strcat(cali_dir,'coeff_exc_exper1.mat'));
load(strcat(data_dir,'\fitting_info_',scan_mode,'.mat'));
fitting_data.coef_det = coeff_det;
fitting_data.coef_exc = coeff_exc;
fitting_data.offset_map = offset;
fitting_data.var_map = var;
fitting_data.gain_map = gain;
[camsize_x,camsize_y] = size(gain);
fitting_info = struct( ...
    'slice_num', 0, ...
    'seg_num', 0, ...% Not more than 4000
    'fit_offset', 0, ...
    'seg_size_xy', size_xy, ...
    'cam_map_size_x', camsize_x, ...
    'cam_map_size_y', camsize_y);
% fit config 1 data
config_intm.slice = 5;
config_intm.fit_os = fit_offset;
config_intm.data = seg_5;
config_intm.seg_config = 1;
config_intm.map_ptr_x = ptr_x_5;
config_intm.map_ptr_y = ptr_y_5;
if fitting_round>1
    fitting_init_config1 = fitting_init(:,seg_5_idx);
    config_intm.fitting_para = reshape(fitting_init_config1,[6*length(ptr_x_5),1]);
else
    config_intm.fitting_para = zeros(6*length(ptr_x_5),1);
end
tic
fitting_res_seg_5 = fit_interface(fitting_info,fitting_data,config_intm);
toc
fitting_results(:,seg_5_idx) = fitting_res_seg_5;
% fit config 2 data
config_intm.slice = 4;
config_intm.fit_os = fit_offset;
config_intm.data = seg_4_down;
config_intm.map_ptr_x = ptr_x_4_down;
config_intm.map_ptr_y = ptr_y_4_down;
if fitting_round>1
    fitting_init_config2 = fitting_init(:,seg_4down_idx);
    config_intm.fitting_para = reshape(fitting_init_config2,[6*length(ptr_x_4_down),1]);
else
    config_intm.fitting_para = zeros(6*length(ptr_x_4_down),1);
end
tic
fitting_res_seg_4down = fit_interface(fitting_info,fitting_data,config_intm);
toc
fitting_results(:,seg_4down_idx) = fitting_res_seg_4down;
% fit config 3 data
config_intm.slice = 4;
config_intm.fit_os = fit_offset;
config_intm.data = seg_4_up;
config_intm.map_ptr_x = ptr_x_4_up;
config_intm.map_ptr_y = ptr_y_4_up;
if fitting_round>1
    fitting_init_config3 = fitting_init(:,seg_4up_idx);
    config_intm.fitting_para = reshape(fitting_init_config3,[6*length(ptr_x_4_up),1]);
else
    config_intm.fitting_para = zeros(6*length(ptr_x_4_up),1);
end
tic
fitting_res_seg_4up = fit_interface(fitting_info,fitting_data,config_intm);
toc
fitting_results(:,seg_4up_idx) = fitting_res_seg_4up;
save(strcat(data_path,'\fitting_result_SM_',scan_mode,'_round',num2str(fitting_round),'.mat'),"fitting_results")
