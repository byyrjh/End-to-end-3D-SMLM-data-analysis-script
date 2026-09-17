clear;clc;close all;
time_binsize=60;
binsize = 1000; % nm
pixel_size_z=10;  % nm
pixel_size_xy=100;  % nm
step_size=400;  % nm
exposuretime=20;  %ms
num_photon_factor=30;
slice_per_stack = 50;
datasave_path='M:\Hao\2022\2022_11_18\cell1\segment_data\stepwise fitting\';
load(strcat(datasave_path,'map_ptr_t_SM_m0.mat'));
load(strcat(datasave_path,'map_ptr_x_SM_m0.mat'));
load(strcat(datasave_path,'map_ptr_y_SM_m0.mat'));
load(strcat(datasave_path,'map_ptr_z_SM_m0.mat'));
load(strcat(datasave_path,'fitting_result_SM_m0.mat'));
load(strcat(datasave_path,'ChiSq_SM_m0.mat'));
load(strcat(datasave_path,'SM_os_aver.mat'));
SM_os_aver = SM_os_aver*pixel_size_z;
SM_os_aver = reshape(SM_os_aver, [27,27,3800]);
x_pos_SM=double((map_ptr_x_SM+fitting_results(1,:)')*pixel_size_xy);
y_pos_SM=double((map_ptr_y_SM+fitting_results(2,:)')*pixel_size_xy);
z_pos_SM=double(map_ptr_z_SM*step_size+fitting_results(3,:)'*pixel_size_z);
LS_os = double(fitting_results(6,:)'*pixel_size_z);
inten_SM=double(fitting_results(4,:)'*num_photon_factor);
t_SM=double(map_ptr_t_SM);
nan_x=isnan(x_pos_SM);
nan_y=isnan(y_pos_SM);
nan_z=isnan(z_pos_SM);
y_oulier = (y_pos_SM<0);
x_oulier = (x_pos_SM<0);
z_oulier = (z_pos_SM<0);
kick_out_idx=nan_x|nan_y|nan_z|y_oulier|x_oulier|z_oulier;%|crlb_idx;
x_pos_SM(kick_out_idx)=[];
y_pos_SM(kick_out_idx)=[];
z_pos_SM(kick_out_idx)=[];
inten_SM(kick_out_idx)=[];
LS_os(kick_out_idx)=[];
t_SM(kick_out_idx)=[];
x_boundry = ceil(max(x_pos_SM));
y_boundry = ceil(max(y_pos_SM));
x_size = ceil(x_boundry/binsize);
y_size = ceil(y_boundry/binsize);
z_size = max(t_SM);
mat_os = zeros(x_size,y_size,z_size);
mat_os_count = zeros(x_size,y_size,z_size);
vol_exposure_t=exposuretime*slice_per_stack*2/1000;
duration_meas=vol_exposure_t*z_size/60; % min
t_axis=linspace(0,duration_meas-exposuretime*slice_per_stack*2/1000/60,z_size);
for i=1:z_size
    idx=(t_SM==i);
    cur_x = x_pos_SM(idx);
    cur_y = y_pos_SM(idx);
    cur_os = LS_os(idx);
    for j = 1:size(cur_x)
        x_cord = ceil(cur_x(j)/binsize);
        y_cord = ceil(cur_y(j)/binsize);
        mat_os(x_cord,y_cord,i) = mat_os(x_cord,y_cord,i)+cur_os(j);
        mat_os_count(x_cord,y_cord,i) = mat_os_count(x_cord,y_cord,i)+1;
    end
end
idx = find(mat_os_count==0);
mat_os_count(idx) = mat_os_count(idx)+1;
mat_os = mat_os./mat_os_count;
map_os = zeros(x_size,y_size);
map_os_ccode = zeros(x_size,y_size);
for i = 1:x_size
    for j = 1:y_size
        LS_os_t_tr = reshape(mat_os(i,j,:),[z_size,1]);
        LS_os_t_tr_ccode = reshape(SM_os_aver(i,j,:),[z_size,1]);
        LS_os_t_tr = tr_smo(LS_os_t_tr,time_binsize);
        LS_os_t_tr_mean = LS_os_t_tr;
        LS_os_t_tr_ccode(isnan(LS_os_t_tr_ccode)) = [];
        map_os_ccode(i,j) = mean(LS_os_t_tr_ccode);
        LS_os_t_tr_mean(isnan(LS_os_t_tr_mean)) = [];
        map_os(i,j) = mean(LS_os_t_tr_mean);
        mat_os(i,j,:) = reshape(LS_os_t_tr,[1,1,z_size]);
    end
end

clim = [-800 800];
figure
imagesc(map_os,clim)
axis equal;
xlim([1 x_size])
ylim([1 y_size])
figure
imagesc(map_os_ccode,clim)
axis equal;
xlim([1 x_size])
ylim([1 y_size])


inpect_x = 17;
inpect_y = 13;
Fs = 1;
fWavelet = 2;
bWidth = 1;
full_x = 1:1:z_size;
test = reshape(mat_os(inpect_x,inpect_y,:),[(z_size),1]);
test_c = reshape(SM_os_aver(inpect_x,inpect_y,:),[(z_size),1]);
sample_x = t_axis(~isnan(test));
test = test(~isnan(test));
sample_x_c = t_axis(~isnan(test_c));
test_c = test_c(~isnan(test_c));
test_interp = interp1(sample_x,test,full_x);
% test_interp_filt = wFilter(test_interp,Fs,fWavelet,bWidth,0);
figure 
hold on
plot(sample_x,test,'LineWidth',1)
% plot(full_x,smooth(test_interp,30))
plot(sample_x_c,test_c,'LineWidth',1)
%ylim([-800 800])
legend("matlab implement","c implement")
box on 
set(gca,'FontSize',18)
xlabel('Time / min')
ylabel('LS offset / nm')


function smoothed_tr = tr_smo(data_tr,binsize)
length_tr = length(data_tr);
for i = 1:length_tr
    if (i-binsize/2<1)
        tr_seg = data_tr(1:i+binsize/2);
    elseif (i+binsize/2>length_tr)
        tr_seg = data_tr(i-binsize/2:length_tr);
    else
        tr_seg = data_tr(i-binsize/2:i+binsize/2);
    end
    tr_seg(tr_seg==0) = [];
    smoothed_tr(i) = mean(tr_seg);
end
end
