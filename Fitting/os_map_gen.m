clear;clc;close all;
scan_mode = 'm0';
scriptFolder = fileparts(mfilename('fullpath'));
packageFolder = fileparts(scriptFolder);
data_dir = fullfile(packageFolder,'data\segment_data\');
cali_dir = fullfile(packageFolder,'data\setup_calibration\');
load(strcat(cali_dir,'camera1_cali.mat'));
load(strcat(data_dir,'fitting_result_SM_',scan_mode,'_round1.mat'));
load(strcat(data_dir,'map_ptr_t_SM_',scan_mode,'.mat'));
load(strcat(data_dir,'map_ptr_z_SM_',scan_mode,'.mat'));
load(strcat(data_dir,'map_ptr_x_SM_',scan_mode,'.mat'));
load(strcat(data_dir,'map_ptr_y_SM_',scan_mode,'.mat'));
load(strcat(data_dir,'seg_config_RMSE_optimal_SM_',scan_mode,'.mat'));
% Center position is relative to the 3rd slice regardless of number of
% slices per segment 
% Therefore granular z position in case 2 has to be corrected by adding 1 slice
bins = 10;
idx_case2 = find(seg_config_RMSE==2);
map_ptr_z_SM_cor = map_ptr_z_SM;
map_ptr_z_SM_cor(idx_case2) = map_ptr_z_SM_cor(idx_case2)+1;
real_z = round(map_ptr_z_SM_cor + fitting_results(3,:)'*10/400);
real_x = ceil((map_ptr_x_SM+fitting_results(1,:)')/bins);
real_y = ceil((map_ptr_y_SM+fitting_results(2,:)')/bins);
bin_edge = -200:5:200;
figure
hold on
h = histogram(fitting_results(6,:),bin_edge);
%% light sheet offset fitting quality control (good data screen)
bin_value = h.Values;
bin_edge = h.BinEdges;
L = tril(ones(length(bin_value),length(bin_value)));
accum1 = L * bin_value';
accum2 = L' * bin_value';
per1 = accum1/max(accum1);
per2 = accum2/max(accum2);
p_cutoff = 0.1;
per1_cut = find(abs(per1-p_cutoff)==min(abs(per1-p_cutoff)));
per2_cut = find(abs(per2-p_cutoff)==min(abs(per2-p_cutoff)));
cut1 = bin_edge(per1_cut);
cut2 = bin_edge(per2_cut);
bins = 10;
[size_xy ~] = size(gain);
size_xy_4D = ceil(size_xy/bins);
size_z_4D = 50;
size_t_4D = max(map_ptr_t_SM);
mat_4D = zeros(size_xy_4D,size_xy_4D,size_z_4D,size_t_4D);
counter_4D = zeros(size_xy_4D,size_xy_4D,size_z_4D,size_t_4D);
idx1 = fitting_results(6,:)>cut1;
idx2 = fitting_results(6,:)<cut2;
idx = idx1&idx2;
os = fitting_results(6,idx);
real_z = real_z(idx);
real_z = min(real_z,50);
real_z = max(1,real_z);
real_x = real_x(idx);
real_y = real_y(idx);
real_t = map_ptr_t_SM(idx);
%% 4D offset map generating
for i = 1:length(real_z)
    mat_4D(real_x(i),real_y(i),real_z(i),real_t(i)) = mat_4D(real_x(i),real_y(i),real_z(i),real_t(i))+os(i);
    counter_4D(real_x(i),real_y(i),real_z(i),real_t(i)) = counter_4D(real_x(i),real_y(i),real_z(i),real_t(i))+1;
end
counter_4D = max(counter_4D,ones(size_xy_4D,size_xy_4D,size_z_4D,size_t_4D));
mat_4D = mat_4D./counter_4D;
%% smooth in t
slide_w = 100;
slide_filter = ones(slide_w,1);
for i = 1:size_xy_4D
    for j = 1:size_xy_4D
        for k = 1:size_z_4D
            temp_data = squeeze(mat_4D(i,j,k,:));
            counter = double(temp_data ~= 0);
            temp_data = conv(temp_data,slide_filter);
            counter = conv(counter,slide_filter);
            counter = max(counter,1);
            temp_data = temp_data./counter;
            mat_4D(i,j,k,:) = temp_data(slide_w/2:end-slide_w/2);
        end
    end
end

%% smoothing in 3D
sigma_f = 2;
dist_ele1 = sqrt(1 + 0.4^2);
dist_ele2 = sqrt(2 + 0.4^2);
dist_ele3 = 0.4;
x = -2:1:2;
[xx yy] = meshgrid(x,x);
[phi rho] = cart2pol(xx,yy);
% filter_mask_1 = [sqrt(2) 1 sqrt(2);
%                1       0      1 ;
%                sqrt(2) 1 sqrt(2)];
% filter_mask_2 = [dist_ele2 dist_ele1 dist_ele2;
%                  dist_ele1 dist_ele3 dist_ele1 ;
%                  dist_ele2 dist_ele1 dist_ele2];
filter_mask_1 = rho;
filter_mask_2 = sqrt(rho.^2 + 0.4^2);
filter_mask_3D(:,:,1) = filter_mask_2;
filter_mask_3D(:,:,2) = filter_mask_1;
filter_mask_3D(:,:,3) = filter_mask_2;
filter_mask_3D = exp(-filter_mask_3D.^2/2/sigma_f^2);
mat_4D_smo = zeros(size(mat_4D));
for i = 1:size_t_4D
    data_3D = mat_4D(:,:,:,i);
    counter_3D = double(data_3D ~= 0);
    data_3D_conv = convn(data_3D,filter_mask_3D);
    counter_3D_conv = convn(counter_3D,filter_mask_3D);
    data_3D_conv = data_3D_conv./counter_3D_conv;
    mat_4D_smo(:,:,:,i) = data_3D_conv(3:end-2,3:end-2,2:end-1);
end
% [plane_fit_para, plane_fit_map_res] = ls_plane_fit(mat_4D_smo);
% save(strcat('plane_fit_para_m',scan_dir,'.mat'),'plane_fit_para');
LS_os_AO_4D_map = reshape(mat_4D_smo,[size_xy_4D*size_xy_4D*size_z_4D,size_t_4D]);
% save(strcat(data_dir,'plane_fit_map_res_',scan_mode,'.mat'),'plane_fit_map_res');
save(strcat(data_dir,'LS_os_AO_4D_map_',scan_mode,'.mat'),'LS_os_AO_4D_map');
%% offset reassign
pos_z = round(map_ptr_z_SM_cor + fitting_results(3,:)'*10/400);
pos_x = ceil((map_ptr_x_SM+fitting_results(1,:)')/bins);
pos_y = ceil((map_ptr_y_SM+fitting_results(2,:)')/bins);
pos_t = map_ptr_t_SM;
fitting_init = fitting_results;
for i = 1:length(map_ptr_t_SM)
    try
        fitting_init(6,i) = mat_4D_smo(pos_x(i),pos_y(i),pos_z(i),pos_t(i));
    catch
    end
end
figure
histogram(fitting_init(6,:),bin_edge)
save(strcat(data_dir,'fitting_init_SM_',scan_mode,'.mat'),'fitting_init');