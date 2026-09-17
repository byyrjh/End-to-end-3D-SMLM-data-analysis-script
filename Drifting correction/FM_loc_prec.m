clear;clc;close all;
time_binsize=30;
FM_outlier = [2 5 10 13 14 15 17 18];% 
piecepos = 1:1:59;% unit hyperstack
flag_savefigure = false;
smooth_f=6;% localization precision is independent of smooth binning size  % 30  % 6
exposuretime=20;  %ms %  
slice_per_stack=50; 
scan_dir = 'm0';  
datasave_path='N:\Lucas\10_LLSM_CudaFitter\03_HeLa_SMLM\2023-08-16-Hela-aTUB\2023-08-16_15-41-01_Hela-ATUB-230815-AB-FOV7-002\segment_data_protect_on_2025_10_1\';

outlier_thres_xy=120;
outlier_thres_z=120;
pixel_size_z=10;  % nm
pixel_size_xy=100;  % nm
step_size=400;  % nm
xybinsize_SM = 1000;  % unit nm
num_photon_factor=30;
% color_map_list = [1 0 0;0 1 0;0 0 1;0 1 1;1 0 1;1 1 0;0 0 0;
%     0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;
%     0, 0, 0.5;0, 0.5, 0;0.5, 0, 0;0, 0.5, 0.5;0.5, 0, 0.5;0.5, 0.5, 0;0.5, 0.5, 0.5
%     0, 0, 0.75;0, 0.75, 0;0.75, 0, 0;0, 0.75, 0.75;0.75, 0, 0.75;0.75, 0.75, 0;0.75, 0.75, 0.75];
mkdir(strcat(datasave_path,scan_dir,'_Traj_analysis'))
datasave_path_traj_ana = strcat(datasave_path,scan_dir,'_Traj_analysis\');
load(strcat(datasave_path,'fitting_info_',scan_dir,'.mat'));
% stationary_pos = stationary_pos-1;
load(strcat(datasave_path,'map_ptr_t_FM_',scan_dir,'.mat'));% M:\Hao\2022\2022_12_11\Fixed HeLa 266nM CLC-Alexa 647 immuno labeled\cell3\segment_data
load(strcat(datasave_path,'map_ptr_t_SM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_x_FM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_y_FM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_z_FM_',scan_dir,'.mat'));
load(strcat(datasave_path,'fitting_result_FM_',scan_dir,'.mat'));
load(strcat(datasave_path,'ChiSq_FM_',scan_dir,'.mat'));
% load(strcat(datasave_path,'SM_os_aver_',scan_dir,'.mat'))
%% fiducial marker track
trace_x_m0=[];
trace_y_m0=[];
trace_z_m0=[];
trace_inten_m0=[];

x_pos=(map_ptr_x_FM+fitting_results(1,:)')*pixel_size_xy;
y_pos=(map_ptr_y_FM+fitting_results(2,:)')*pixel_size_xy;
% z_pos=map_ptr_z_FM*step_size+fitting_results(3,:)'*pixel_size_z;
z_pos=map_ptr_z_FM*step_size+fitting_results(3,:)'*pixel_size_z;
ls_os_in_fit = fitting_results(6,:)'*pixel_size_z;
fiducial_idx = 1:1:length(x_pos);
cor_lab=[x_pos y_pos z_pos];

data_full=[cor_lab fitting_results(4,:)'*num_photon_factor ls_os_in_fit fitting_results(3,:)' fiducial_idx' map_ptr_t_FM];
data_full=double((data_full));

maxdisp=1000;  % nm
param=struct('mem',0,'dim',3,'good',1,'quiet',1);

track_output_m0=track(data_full,maxdisp,param);
% data structure
% x y z brightness offset z_fitting fiducial_index t trace_index
% 1-5 entries general param
% 6 to investigate digitization
% 7 to compare lightsheet offset fitting by 1D Gaussian fitting
traj_m0=zeros(max(track_output_m0(:,9)),1);

for i=1:size(track_output_m0,1)
    traj_m0(track_output_m0(i,9))=traj_m0(track_output_m0(i,9))+1;
end

marker_idx_m0=find(traj_m0==max(map_ptr_t_FM));
origin_arr=zeros(length(marker_idx_m0),3);
%% raw trajactory generation
trace_x_m0=[];
trace_y_m0=[];
trace_z_m0=[];
trace_inten_m0=[];
trace_os_m0 = [];
trace_z_seg = [];
trace_idx = [];
for i=1:length(marker_idx_m0)
    temp_x=track_output_m0(track_output_m0(:,9)==marker_idx_m0(i),1);
    temp_y=track_output_m0(track_output_m0(:,9)==marker_idx_m0(i),2);
    temp_z=track_output_m0(track_output_m0(:,9)==marker_idx_m0(i),3);
    temp_inten=track_output_m0(track_output_m0(:,9)==marker_idx_m0(i),4);
    temp_os = track_output_m0(track_output_m0(:,9)==marker_idx_m0(i),5);
    temp_z_seg = track_output_m0(track_output_m0(:,9)==marker_idx_m0(i),6);
    temp_idx = track_output_m0(track_output_m0(:,9)==marker_idx_m0(i),7);
    trace_x_m0=[trace_x_m0 temp_x];
    trace_y_m0=[trace_y_m0 temp_y];
    trace_z_m0=[trace_z_m0 temp_z];
    origin_arr(i,1)=mean(temp_x(1:30));
    origin_arr(i,2)=mean(temp_y(1:30));
    origin_arr(i,3)=mean(temp_z(1:30));
    trace_inten_m0=[trace_inten_m0 temp_inten];
    trace_os_m0 = [trace_os_m0 temp_os];
    trace_z_seg = [trace_z_seg temp_z_seg];
    trace_idx = [trace_idx temp_idx];
end
[vol_num,~]=size(trace_inten_m0);
vol_idx=linspace(1,vol_num,vol_num);
marker_idx_m0(FM_outlier) = [];
trace_x_m0(:,FM_outlier) = [];
trace_y_m0(:,FM_outlier) = [];
trace_z_m0(:,FM_outlier) = [];
trace_inten_m0(:,FM_outlier) = [];
origin_arr(FM_outlier,:) = [];
trace_os_m0(:,FM_outlier) = [];
trace_z_seg(:,FM_outlier) = [];
trace_idx(:,FM_outlier) = [];
num_vol=length(trace_inten_m0);
vol_exposure_t=exposuretime*slice_per_stack*2/1000;
duration_meas=vol_exposure_t*num_vol/60; % min
t_axis=linspace(0,duration_meas-exposuretime*slice_per_stack*2/1000/60,num_vol);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Traces are sorted so z in descending order %%%
[~, z_rank_plot] = sort(origin_arr(:,3)'/1000);
trace_x_m0 = trace_x_m0(:,z_rank_plot);
trace_y_m0 = trace_y_m0(:,z_rank_plot);
trace_z_m0 = trace_z_m0(:,z_rank_plot);
trace_inten_m0 = trace_inten_m0(:,z_rank_plot);
origin_arr = origin_arr(z_rank_plot,:);
trace_os_m0 = trace_os_m0(:,z_rank_plot);
trace_z_seg = trace_z_seg(:,z_rank_plot)*pixel_size_z;
trace_idx = trace_idx(:,z_rank_plot);
color_map_list = cm_gen(length(marker_idx_m0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This block is used to compared lightsheet offset fitting by fitter and 1D Gaussian fitting
raw_z_smooth = smooth_gen(trace_z_seg, vol_per_hyper, piecepos, 20);
figure
hold on
for i = 1:size(trace_z_seg,2)
    plot(t_axis,raw_z_smooth(:,i),'Color',color_map_list(i,:))
end
figure 
hold on 
for i = 1:size(trace_z_seg,2)
    plot(t_axis,trace_os_m0(:,i),'Color',color_map_list(i,:))
end

os_to_cen = raw_z_smooth + trace_os_m0;
[temp_x temp_y] = size(os_to_cen);
os_vec = reshape(os_to_cen,[temp_x*temp_y 1]);
idx_vec = reshape(trace_idx,[temp_x*temp_y 1]);
load(strcat(datasave_path,'mu_vec.mat'));
lsfit_Gau = mu_vec(idx_vec);
% err_lsos_fit = lsfit_Gau - os_vec;
lsfit_Gau_mat = reshape(lsfit_Gau,[temp_x temp_y]);
lsfit_Gau_smooth = smooth_gen(lsfit_Gau_mat, vol_per_hyper, piecepos, 20);
err_lsos_fit = reshape(lsfit_Gau_smooth,[temp_x*temp_y 1]) - os_vec;
errrrrr_lsos_fit = lsfit_Gau - os_vec;
bin_edge = -400:5:400;
figure
hold on
histogram(errrrrr_lsos_fit,bin_edge)
histogram(err_lsos_fit,bin_edge)

%% PLOT
% x y z raw
figure
hold on
xlabel('Time / min')
ylabel('Drifting / nm')
box on 
grid on
title('Drifting in x (raw)')
set(gca,'FontSize',20)
mean_trace_x = [];

for i=1:length(marker_idx_m0)
    plot(t_axis,trace_x_m0(:,i)-trace_x_m0(1,i),'Color',color_map_list(i,:))
    mean_trace_x = [mean_trace_x trace_x_m0(:,i)-trace_x_m0(1,i)];
end
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Drifting in x_raw.bmp'))
end

figure
hold on
xlabel('Time / min')
ylabel('Drifting / nm')
box on 
grid on
title('Drifting in propagation direction (raw)')
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_y_m0(:,i)-trace_y_m0(1,i),'Color',color_map_list(i,:))
end
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Drifting in y_raw.bmp'))
end

figure
hold on
xlabel('Time / min')
ylabel('Drifting / nm')
box on 
grid on
title('Drifting in z (raw)')
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_z_m0(:,i)-trace_z_m0(1,i),'Color',color_map_list(i,:))
end
trace_z_abs = trace_z_m0 - trace_z_m0(1,:);
% digitized localization (trace_z_seg) is possibly due to local
% nonmonotonic feature of calibrated PSF in z 
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Drifting in z_raw.bmp'))
end
%% PLOT 
% photobleaching
figure
hold on
title('Photobleaching')
xlabel('Time / min')
ylabel('Number of photons')
box on 
grid on
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_inten_m0(:,i),'Color',color_map_list(i,:))
end
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Photobleaching.bmp'))
end
%% PLOT
% light sheet drifting trajactory
% vol_per_hyper = vol_per_hyper-stationary_pos;
num_hyperstack = num_vol/vol_per_hyper;
tot_vol = vol_per_hyper*num_hyperstack;
vol_time = exposuretime*slice_per_stack*2;
hyper_pos = linspace(0,tot_vol*vol_time-vol_time*vol_per_hyper,num_hyperstack)/1000/60;

load(strcat(datasave_path,'LS_os_',scan_dir,'.mat'));
figure
hold on
title('Light sheet offset FM')
xlabel('Time / min')
ylabel('$\bar{\textrm{z}}_\textrm{mf}$ / nm','Interpreter','latex')
box on 
grid on
set(gca,'FontSize',12)

for i=1:length(marker_idx_m0)
    plot(t_axis,trace_os_m0(:,i),'LineWidth',1.5,'Color',color_map_list(i,:));
end
for i=1:num_hyperstack
    if mod(i,10) == 0
        xline(hyper_pos(i),'LineWidth',1)
    else
        xline(hyper_pos(i),'--','LineWidth',0.5)
    end
end
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Light sheet offset.bmp'))
end

trace_os_SM = [];
x_range = ceil(camsize_x / (xybinsize_SM / pixel_size_xy));
for i=1:length(marker_idx_m0)
    temp_x = ceil(origin_arr(i,1)/pixel_size_xy/(xybinsize_SM / pixel_size_xy));
    temp_y = ceil(origin_arr(i,2)/pixel_size_xy/(xybinsize_SM / pixel_size_xy));
    idx_bin_pix = (temp_y-1)*x_range+temp_x;
%     trace_os_SM = [trace_os_SM;SM_os_aver(idx_bin_pix,:)];
end
%% smoothed trajectory generation
drifting_corr_x = smooth_gen(trace_x_m0, vol_per_hyper, piecepos, smooth_f);
drifting_corr_y = smooth_gen(trace_y_m0, vol_per_hyper, piecepos, smooth_f);
drifting_corr_z = smooth_gen(trace_z_m0, vol_per_hyper, piecepos, smooth_f);
% 7.2.2025 digitalization of z localization originates from over fine z
% step size of calibrated PSF_det (10 nm), which is not monotonically
% within the optimization domain that cases ringing in downhill converging
% cost function

drift_smooth.x = [];
drift_smooth.y = [];
drift_smooth.z = [];


%% PLOT
% smoothed x y 
figure
hold on
grid on
box on
set(gca,'FontSize',20)
title('Drifting in x (mean)')
for i=1:length(marker_idx_m0)
    plot(t_axis,drifting_corr_x(:,i)-drifting_corr_x(1,i),'Color',color_map_list(i,:),'LineWidth',1)
    drift_smooth.x = [drift_smooth.x drifting_corr_x(:,i)-drifting_corr_x(1,i)];
end
for i=1:num_hyperstack
    if mod(i,10) == 0
        xline(hyper_pos(i),'LineWidth',1)
    else
        xline(hyper_pos(i),'--','LineWidth',0.5)
    end
end
xlabel('Time / min')
ylabel('Drifting / nm')
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Drifting in x_mean.bmp'))
end

figure
hold on
grid on
box on
set(gca,'FontSize',20)
title('Drifting in propagation direction')
for i=1:length(marker_idx_m0)
    plot(t_axis,drifting_corr_y(:,i)-drifting_corr_y(1,i),'Color',color_map_list(i,:),'LineWidth',1);
    drift_smooth.y = [drift_smooth.y drifting_corr_y(:,i)-drifting_corr_y(1,i)];
end
for i=1:num_hyperstack
    if mod(i,10) == 0
        xline(hyper_pos(i),'LineWidth',1)
    else
        xline(hyper_pos(i),'--','LineWidth',0.5)
    end
end
xlabel('Time / min')
ylabel('Drifting / nm')
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Drifting in propagation_mean.bmp'))
end
%% PLOT
% z trace without no correction
figure
hold on
grid on
box on
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    % lsos_smooth = trace_os_m0(:,i)-trace_os_m0(1,i);
    plot(t_axis,drifting_corr_z(:,i)-drifting_corr_z(1,i),'Color',color_map_list(i,:),'LineWidth',1)
    % drift_smooth.z = [drift_smooth.z drifting_corr_z(:,i)-drifting_corr_z(1,i)+lsos_smooth];
end
drifting_err_z_wo_cor = drifting_corr_z - drifting_corr_z(1,:);
drifting_err_z_wo_cor = std(drifting_err_z_wo_cor,0,2);
xlabel('Time / min')
ylabel('$\mathrm{z}_{\mathrm{mf}} $ / nm','Interpreter','latex')
% title('Drifting in z (mean) wo cor')

for i=1:num_hyperstack
    if mod(i,10) == 0
        xline(hyper_pos(i),'LineWidth',1)
    else
        xline(hyper_pos(i),'--','LineWidth',0.5)
    end
end

if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Drifting in z (mean) wo cor.bmp'))
end

%% PLOT
% plot focal plane corrected z (intermediate)
figure
hold on
grid on
set(gca,'FontSize',20)
drifting_corr_z_temp = drifting_corr_z + trace_os_m0;
drifting_corr_z_temp = drifting_corr_z_temp - drifting_corr_z_temp(1,:);
drifting_err_z_fp_cor = std(drifting_corr_z_temp,0,2);
for i=1:length(marker_idx_m0)
    % lsos_smooth = trace_os_m0(:,i)-trace_os_m0(1,i);
    % drift_smooth_temp = drifting_corr_z(:,i) + trace_os_m0(:,i);
%     plot(t_axis,drifting_corr_z(:,i)-drifting_corr_z(1,i)+lsos_smooth,'Color',color_map_list(i,:),'LineWidth',1)
    plot(t_axis,drifting_corr_z_temp(:,i),'Color',color_map_list(i,:),'LineWidth',1)
    % drift_smooth.z = [drift_smooth.z drifting_corr_z(:,i)-drifting_corr_z(1,i)+lsos_smooth];
end

for i=1:num_hyperstack
    if mod(i,10) == 0
        xline(hyper_pos(i),'LineWidth',1)
    else
        xline(hyper_pos(i),'--','LineWidth',0.5)
    end
end
xlabel('Time / min')
ylabel('$\mathrm{z}_{\mathrm{mf}} + \bar{\textrm{z}}_\textrm{mf} $ / nm','Interpreter','latex')
% title('Drifting in z (mean) w FP cor') %focal plane
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Drifting in z (mean) w FP cor.bmp'))
end

%% light sheet plane fitting and plotting
mat_coeff = zeros(vol_num,3);
% model ax+by+c = z
for i=1:vol_num
    LS_plane_data = trace_os_m0(i,:);
    x_bin_pos = drifting_corr_x(i,:)/100/10; % 100 is pixel size 10 is bin size
    y_bin_pos = drifting_corr_y(i,:)/100/10; % 100 is pixel size 10 is bin size
    A = [x_bin_pos-1; y_bin_pos-1; ones(size(y_bin_pos))]';
    temp_coeff = inv(A'*A)*(A')*(LS_plane_data');
    mat_coeff(i,:) = temp_coeff';
end
for i = 1:3
    mat_coeff(:,i) = smooth(mat_coeff(:,i),50);
end
f = figure;
subplot(1,3,1)
plot(t_axis,mat_coeff(:,1))
xlabel('Time / min')
ylabel('scan slope / nm/m')
set(gca,'FontSize',20)
subplot(1,3,2)
plot(t_axis,mat_coeff(:,2))
xlabel('Time / min')
ylabel('propagation slope / nm/m')
set(gca,'FontSize',20)
subplot(1,3,3)
plot(t_axis,mat_coeff(:,3))
xlabel('Time / min')
ylabel('Intercept / nm')
set(f, 'Position',[100 100 1200 400]);
set(gca,'FontSize',20)

if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'light sheet plane fitting.bmp'))
end

trace_os_m0_fitted = [];
for i=1:size(drifting_corr_x,2)
    x_bin_pos = drifting_corr_x(:,i)/100/10-1; % 100 is pixel size 10 is bin size
    y_bin_pos = drifting_corr_y(:,i)/100/10-1; % 100 is pixel size 10 is bin size
    z_fitted = x_bin_pos.*mat_coeff(:,1) + y_bin_pos.*mat_coeff(:,2) + mat_coeff(:,3);
    trace_os_m0_fitted = [trace_os_m0_fitted z_fitted];
end
%% correct z with light sheet plane fitting and focal plane drifting
%%
figure
hold on
grid on
set(gca,'FontSize',20)
% rel_lsos = trace_os_m0_fitted - trace_os_m0_fitted(1,:); % relative light sheet mechanical drifting over time
% lsos_fp_fluc = trace_os_m0 - rel_lsos;
% lsos_fp_fluc = lsos_fp_fluc - lsos_fp_fluc(1,:); 
drifting_err_z_fp_cor_ls_mech_cor_trace = drifting_corr_z + trace_os_m0 - trace_os_m0_fitted;
drifting_err_z_fp_cor_ls_mech_cor_trace = drifting_err_z_fp_cor_ls_mech_cor_trace - drifting_err_z_fp_cor_ls_mech_cor_trace(1,:);
% drifting_err_z_fp_cor_ls_mech_cor_temp = drifting_err_z_fp_cor_ls_mech_cor_temp - ...
%     drifting_err_z_fp_cor_ls_mech_cor_temp(1,:);
drifting_err_z_fp_cor_ls_mech_cor = std(drifting_err_z_fp_cor_ls_mech_cor_trace,0,2);
for i=1:length(marker_idx_m0)
    % plot(t_axis,drifting_err_z_fp_cor_ls_mech_cor_temp(:,i),'Color',color_map_list(i,:),'LineWidth',1)
     plot(t_axis,drifting_err_z_fp_cor_ls_mech_cor_trace(:,i),'Color',color_map_list(i,:),'LineWidth',1)
end

for i=1:num_hyperstack
    if mod(i,10) == 0
        xline(hyper_pos(i),'LineWidth',1)
    else
        xline(hyper_pos(i),'--','LineWidth',0.5)
    end
end
xlabel('Time / min')
ylabel('$\mathrm{z}_{\mathrm{mf}} + \bar{\textrm{z}}_\textrm{mf} - \bar{\textrm{z}}_\textrm{mf}^\textrm{fit} $ / nm','Interpreter','latex')
% title('Drifting in z (mean) w FP cor with ls drift cor') %focal plane
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Drifting in z (mean) w FP cor with ls drift cor.bmp'))
end

%% No idea what was my thought back then
%{
[z_cor, z_rank] = sort(origin_arr(:,3)'/1000);
drifting_corr_z_rank = drifting_corr_z(:,z_rank);
drifting_corr_z_rank = drifting_corr_z_rank-drifting_corr_z_rank(1,:);
trace_os_m0_rank = trace_os_m0(:,z_rank);
drifting_corr_z_rank_cor = drifting_corr_z_rank+(trace_os_m0_rank-trace_os_m0_rank(1,:));

rmse_wo_cor = zeros(size(drifting_corr_z_rank(:,1)));
rmse_w_cor = zeros(size(drifting_corr_z_rank(:,1)));
for i = 1:length(drifting_corr_z_rank(:,1))
    coef_wo_cor = polyfit(z_cor,drifting_corr_z_rank(i,:),1);
    coef_w_cor = polyfit(z_cor,drifting_corr_z_rank_cor(i,:),1);
    fit_data_wo_cor = z_cor*coef_wo_cor(1)+coef_wo_cor(2);
    fit_data_w_cor = z_cor*coef_w_cor(1)+coef_w_cor(2);
    rmse_wo_cor(i) = sqrt(mean((fit_data_wo_cor-drifting_corr_z_rank(i,:)).^2));
    rmse_w_cor(i) = sqrt(mean((fit_data_w_cor-drifting_corr_z_rank_cor(i,:)).^2));
end
figure
hold on 
plot(t_axis,rmse_wo_cor,'LineWidth',2)
plot(t_axis,rmse_w_cor,'LineWidth',2)
legend('wo correction','w correction')
box on
grid on
set(gca,'FontSize',18)
xlabel('Time / min')
ylabel('RMSE / nm')
xlim([1 round(max(t_axis))])
%}
%% Comparison between different correction methods
figure
hold on
plot(t_axis,smooth(drifting_err_z_wo_cor,30),'LineWidth',2);
plot(t_axis,smooth(drifting_err_z_fp_cor,30),'LineWidth',2);
plot(t_axis,smooth(drifting_err_z_fp_cor_ls_mech_cor,30),'LineWidth',2);
legend("wo cor","fp cor","fp cor + ls mech cor")
box on 
grid on 
xlabel('Time / min')
ylabel('RMSE / nm')
set(gca,'FontSize',18)
xlim([1 round(max(t_axis))])
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'comp drifting err in z.bmp'))
end
%% PLOT
% fiducial marker scattering
origin_arr = origin_arr/1000;
figure
hold on
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
% scatter3(origin_arr(i,1),origin_arr(i,2),origin_arr(i,3),'LineWidth',4,'MarkerEdgeColor',color_map_list(i,:));
scatter(origin_arr(i,2),origin_arr(i,3),'LineWidth',4,'MarkerEdgeColor',color_map_list(i,:));
end
axis equal
% xlabel('scanning / m')
xlabel('Propagation dir. / m')
ylabel('Detection dir. / m')
box on
title('Fiducial marker distribution')
% xlim([5 20])
% ylim([0 20])
% view(90,0);
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Fiducial marker distribution.bmp'));
end

%% localization precision generation
trace_x_corrected=trace_x_m0-drifting_corr_x;
trace_y_corrected=trace_y_m0-drifting_corr_y;
trace_z_corrected=trace_z_m0-drifting_corr_z;
% kick out spikes
for i=1:length(marker_idx_m0)
    idx_temp=(abs(trace_x_corrected(:,i))>outlier_thres_xy)|(abs(trace_y_corrected(:,i))>outlier_thres_xy)|(abs(trace_z_corrected(:,i))>outlier_thres_z);
    idx_pos=find(idx_temp==true);
    if ~isempty(idx_pos)
        for j=1:length(idx_pos)
            if (idx_pos(j)>smooth_f/2 && (idx_pos(j)+smooth_f/2)<num_vol)
                temp_x=trace_x_m0(idx_pos(j)-smooth_f/2:idx_pos(j)+smooth_f/2,i);
                temp_y=trace_y_m0(idx_pos(j)-smooth_f/2:idx_pos(j)+smooth_f/2,i);
                temp_z=trace_z_m0(idx_pos(j)-smooth_f/2:idx_pos(j)+smooth_f/2,i);
                temp_x(smooth_f/2+1)=[];
                temp_y(smooth_f/2+1)=[];
                temp_z(smooth_f/2+1)=[];
                trace_x_m0(idx_pos(j),i)=mean(temp_x);
                trace_y_m0(idx_pos(j),i)=mean(temp_y);
                trace_z_m0(idx_pos(j),i)=mean(temp_z);
            end
        end
    end
end

drifting_corr_x = smooth_gen(trace_x_m0, vol_per_hyper, piecepos, smooth_f);
drifting_corr_y = smooth_gen(trace_y_m0, vol_per_hyper, piecepos, smooth_f);
drifting_corr_z = smooth_gen(trace_z_m0, vol_per_hyper, piecepos, smooth_f);
trace_x_corrected=trace_x_m0-drifting_corr_x;
trace_y_corrected=trace_y_m0-drifting_corr_y;
trace_z_corrected=trace_z_m0-drifting_corr_z;
%%
figure
hold on
grid on
box on
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_x_corrected(:,i))
end
xlabel('Time / min')
ylabel('Drifting / nm')
title('Corrected localiation traces in x')
ylim([-200 200])
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Corrected localiation traces in x.bmp'))
end

figure
hold on
grid on
box on
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_y_corrected(:,i))
end
xlabel('Time / min')
ylabel('Drifting / nm')
title('Corrected localiation traces in propagation direction')
ylim([-200 200])
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Corrected localiation traces in propagation.bmp'))
end

figure
hold on
grid on
box on
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_z_corrected(:,i))
end
xlabel('Time / min')
ylabel('Drifting / nm')
title('Corrected localiation traces in z')
ylim([-200 200])
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Corrected localiation traces in z.bmp'))
end

loc_precision_x=std(trace_x_corrected,0,1);
loc_precision_y=std(trace_y_corrected,0,1);
loc_precision_z=std(trace_z_corrected,0,1);
aver_loc_x=mean(loc_precision_x);
aver_loc_y=mean(loc_precision_y);
aver_loc_z=mean(loc_precision_z);
std_loc_x=std(loc_precision_x);
std_loc_y=std(loc_precision_y);
std_loc_z=std(loc_precision_z);

aver_loc_x_t=[];
aver_loc_y_t=[];
aver_loc_z_t=[];
for i=1:floor(num_vol/time_binsize)
    aver_loc_x_t=[aver_loc_x_t; (std(trace_x_corrected(1+(i-1)*time_binsize:i*time_binsize,:),0,1))];
    aver_loc_y_t=[aver_loc_y_t; (std(trace_y_corrected(1+(i-1)*time_binsize:i*time_binsize,:),0,1))];
    aver_loc_z_t=[aver_loc_z_t; (std(trace_z_corrected(1+(i-1)*time_binsize:i*time_binsize,:),0,1))];
end
bin_time_axis=(time_binsize*vol_exposure_t:time_binsize*vol_exposure_t:time_binsize*vol_exposure_t*floor(num_vol/time_binsize))/60;
std_x=[];
std_y=[];
std_z=[];
figure
hold on
plot(bin_time_axis,mean(aver_loc_x_t,2),'LineWidth',2)
plot(bin_time_axis,mean(aver_loc_y_t,2),'LineWidth',2)
plot(bin_time_axis,mean(aver_loc_z_t,2),'LineWidth',2)
grid on
box on
set(gca,'FontSize',20)
xlabel('Time / min')
ylabel('Localization precision / nm')
title('Localization precision vs time')
ylim([0 50])
legend('x','y','z')
if flag_savefigure
    saveas(gcf,strcat(datasave_path_traj_ana,'Localization precision vs time.bmp'))
end
%% changing smooth factor
aver_loc_x=[];
aver_loc_y=[];
aver_loc_z=[];
std_loc_x=[];
std_loc_y=[];
std_loc_z=[];
for smooth_factor=2:4:102
    drifting_corr_x = smooth_gen(trace_x_m0, vol_per_hyper, piecepos, smooth_factor);
    drifting_corr_y = smooth_gen(trace_y_m0, vol_per_hyper, piecepos, smooth_factor);
    drifting_corr_z = smooth_gen(trace_z_m0, vol_per_hyper, piecepos, smooth_factor);
    trace_x_corrected = trace_x_m0 - drifting_corr_x;
    trace_y_corrected = trace_y_m0 - drifting_corr_y;
    trace_z_corrected = trace_z_m0 - drifting_corr_z;
    loc_precision_x=std(trace_x_corrected,0,1);
    loc_precision_y=std(trace_y_corrected,0,1);
    loc_precision_z=std(trace_z_corrected,0,1);
    aver_loc_x=[aver_loc_x mean(loc_precision_x)];
    aver_loc_y=[aver_loc_y mean(loc_precision_y)];
    aver_loc_z=[aver_loc_z mean(loc_precision_z)];
    std_loc_x=[std_loc_x std(loc_precision_x)];
    std_loc_y=[std_loc_y std(loc_precision_y)];
    std_loc_z=[std_loc_z std(loc_precision_z)];
end
smooth_factor=2:4:102;
figure
hold on
errorbar(smooth_factor,aver_loc_x,std_loc_x,'LineWidth',2)
errorbar(smooth_factor,aver_loc_y,std_loc_y,'LineWidth',2)
errorbar(smooth_factor,aver_loc_z,std_loc_z,'LineWidth',2)
grid on 
box on
set(gca,'FontSize',20)
xlabel('Sliding window width')
ylabel('Localization precision / nm')
xlim([0 100])
legend('x','y','z')

