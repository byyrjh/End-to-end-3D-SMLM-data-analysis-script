function corrected_data = drifting_correct(data,FM_trace,setting_para,scan)
smooth_f = setting_para.smooth_f;
vol_per_hyper = setting_para.vol_per_hyper-setting_para.stationary_pos;
drifting_corr_x = smooth_gen(FM_trace.x, vol_per_hyper, setting_para.piecepos, smooth_f);
drifting_corr_y = smooth_gen(FM_trace.y, vol_per_hyper, setting_para.piecepos, smooth_f);
drifting_corr_z = smooth_gen(FM_trace.z, vol_per_hyper, setting_para.piecepos, smooth_f) + FM_trace.LS_os;
%%%%%%  light sheet offset fit for fiducial marker  11.11.2025 %%%%%%%%%
drifting_corr_z_mech = ls_mech_drift_FM(FM_trace);
drifting_corr_z = drifting_corr_z - drifting_corr_z_mech;
load(strcat('N:\Lucas\10_LLSM_CudaFitter\03_HeLa_SMLM\2023-08-16-Hela-aTUB\2023-08-16_15-41-01_Hela-ATUB-230815-AB-FOV7-002\segment_data\plane_fit_map_res_',scan,'.mat'));
x_bin = ceil(data.x/1000);
y_bin = ceil(data.y/1000);
z_bin = ceil(data.z/400);
os_res_mech = zeros(size(z_bin));
for i = 1:length(z_bin)
    os_res_mech(i) = plane_fit_map_res(x_bin(i),y_bin(i),z_bin(i),data.t(i));
end
%%%%%%%%%%%%%%%%%%%%%  drifting correction 12.3.2024  %%%%%%%%%%%%%%%%%%%%%
% generate support coordinate z at each t
% fit parabolic curves for x and y
% correct drifting in x and y
% correct drifting in z  
for i = 1:max(data.t)
    support_coordinate_z = drifting_corr_z(i,:);
    sample_x = drifting_corr_x(i,:) - drifting_corr_x(1,:);
    sample_y = drifting_corr_y(i,:) - drifting_corr_y(1,:);
    coef_x = polyfit(support_coordinate_z,sample_x,2);
    coef_y = polyfit(support_coordinate_z,sample_y,2);
    driting_z_mean = mean(drifting_corr_z(i,:)-drifting_corr_z(1,:));
    t_idx=(data.t==i);
    data.x(t_idx)=data.x(t_idx)-(coef_x(1)*data.z(t_idx).^2+coef_x(2)*data.z(t_idx)+coef_x(3));
    data.y(t_idx)=data.y(t_idx)-(coef_y(1)*data.z(t_idx).^2+coef_y(2)*data.z(t_idx)+coef_y(3));
    data.z(t_idx)=data.z(t_idx)-os_res_mech(t_idx)-driting_z_mean;
end
corrected_data = data;
end