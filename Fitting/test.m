clear;clc;
% measured excitation PSF 640 nm, stepsize 10 nm, 0.2 NA sq Lattice
% This script is simplified version of seg_vis.m considering 5 slices fitting
% Compared to seg_vis this script runs faster 
load('psf_exc_wo_bg.mat');
z_exc = linspace(-1200,1200,241);
exc_tot = exc_tot/max(exc_tot);
options = fitoptions('gauss1');
options.Lower = [0.8 -600 745];
options.Upper = [Inf 600 765];
[fit_para_lat1 fg] = fit(z_exc',exc_tot,'gauss1',options);
figure
hold on
plot(z_exc,exc_tot,marker = "*",LineStyle="none");
z_fit = fit_para_lat1.a1*exp(-((z_exc-fit_para_lat1.b1)/fit_para_lat1.c1).^2);
plot(z_exc,z_fit)
sigma_min = fit_para_lat1.c1-10;
sigma_max = fit_para_lat1.c1+10;
load('map_ptr_t_SM_m0.mat');
load('map_ptr_x_SM_m0.mat');
load('map_ptr_y_SM_m0.mat');
load('map_ptr_x_SM_m0.mat');
load('camera1_cali.mat');
load('seg_data_SM_m0.mat');
num_seg = length(map_ptr_x_SM);
inten_tr = zeros(5,num_seg);
for i = 1:num_seg
    temp_seg = seg_data_SM(:,:,(i-1)*5+1:i*5);
    temp_os = offset(map_ptr_x_SM(i)-9:map_ptr_x_SM(i)+9,map_ptr_y_SM(i)-9:map_ptr_y_SM(i)+9);
    temp_gain = gain(map_ptr_x_SM(i)-9:map_ptr_x_SM(i)+9,map_ptr_y_SM(i)-9:map_ptr_y_SM(i)+9);
    temp_seg = (temp_seg-temp_os)./temp_gain;
    inten_tr(:,i) = squeeze(sum(temp_seg,[1 2]));
end

z = -800 : 400 : 800;
z = z';
RMSE_vec = nan(num_seg,1);
Rsq_vec = nan(num_seg,1);
On_config = zeros(num_seg,1);
seg_config_RMSE = zeros(num_seg,1);
seg_config_Rsq = zeros(num_seg,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% case 1  1 1 1 1 1
%%%%%%% case 2  0 1 1 1 1
%%%%%%% case 3  1 1 1 1 0
%%%%%%% case 4  0 0 1 1 1
%%%%%%% case 5  1 1 1 0 0
%%%%%%% case 6  0 1 1 1 0
%%%%%%% case 7  0 0 0 1 1
%%%%%%% case 8  0 0 1 1 0
%%%%%%% case 9  0 1 1 0 0
%%%%%%% case 10 1 1 0 0 0
idx_mat = cell(10,1);
idx_mat{1,1} = logical([1 1 1 1 1]);
idx_mat{2,1} = logical([0 1 1 1 1]);
idx_mat{3,1} = logical([1 1 1 1 0]);
idx_mat{4,1} = logical([0 0 1 1 1]);
idx_mat{5,1} = logical([1 1 1 0 0]);
idx_mat{6,1} = logical([0 1 1 1 0]);
idx_mat{7,1} = logical([0 0 0 1 1]);
idx_mat{8,1} = logical([0 0 1 1 0]);
idx_mat{9,1} = logical([0 1 1 0 0]);
idx_mat{10,1} = logical([1 1 0 0 0]);
z_mat = cell(10,1);
for i = 1:10
    z_mat{i,1} = z(idx_mat{i,1});
end
tic
parfor i = 1 : num_seg
    temp = inten_tr(:,i);
    temp = temp / max(temp);
    x = z_mat{1,1};
    try
        [fit_para_lat1 fg] = fit(x,temp,'gauss1',options);
        RMSE_vec(i) = fg.rmse;
        Rsq_vec(i) = fg.rsquare;
    catch
    end
end
toc
%%%%%%%%%%%%%%%%%
bin_edge_rmse = linspace(0.1,0.6,100);
bin_edge_rsq = linspace(0,1,100);
figure
histogram(RMSE_vec,bin_edge_rmse)
xlabel('Normalized RMSE')
figure
histogram(Rsq_vec,bin_edge_rsq)
xlabel('R-squared')




