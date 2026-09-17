close all;clear;clc;
% measured excitation PSF 640 nm, stepsize 10 nm, 0.2 NA sq Lattice
scriptFolder = fileparts(mfilename('fullpath'));
packageFolder = fileparts(scriptFolder);
data_path = fullfile(packageFolder,'data\segment_data');
load(fullfile(packageFolder,'data\psf_exc_wo_bg.mat'));
z_exc = linspace(-1200,1200,241);
scan_mode = 'm0';
% exc_tot = exc_tot;% No need to normalize
fit_model = 'a1*exp(-((x-b1)/c1)^2)+d1'; % to increase degree of freedom, fix sigma
options = fitoptions(fit_model);
options.Lower = [0.8 -600 745 0];
options.Upper = [Inf 600 765 Inf];
options.StartPoint = [1 0 755 0];
[fit_para_lat1 fg] = fit(z_exc',exc_tot,fit_model,options);
figure
hold on
plot(z_exc,exc_tot,marker = "*",LineStyle="none");
z_fit = fit_para_lat1.a1*exp(-((z_exc-fit_para_lat1.b1)/fit_para_lat1.c1).^2);
plot(z_exc,z_fit)
load(strcat(data_path,'\map_ptr_t_SM_',scan_mode,'.mat'));
load(strcat(data_path,'\map_ptr_x_SM_',scan_mode,'.mat'));
load(strcat(data_path,'\map_ptr_y_SM_',scan_mode,'.mat'));
load(strcat(data_path,'\map_ptr_x_SM_',scan_mode,'.mat'));
load(fullfile(packageFolder,'\data\setup_calibration\camera1_cali.mat'));
load(strcat(data_path,'\seg_data_SM_',scan_mode,'.mat'));
num_seg = length(map_ptr_x_SM);
inten_tr = zeros(5,num_seg);
for i = 1:num_seg
    temp_seg = seg_data_SM(:,:,(i-1)*5+1:i*5);
    temp_os = offset(map_ptr_x_SM(i)-9:map_ptr_x_SM(i)+9,map_ptr_y_SM(i)-9:map_ptr_y_SM(i)+9);
    temp_gain = gain(map_ptr_x_SM(i)-9:map_ptr_x_SM(i)+9,map_ptr_y_SM(i)-9:map_ptr_y_SM(i)+9);
    temp_seg = (temp_seg-temp_os)./temp_gain;
    inten_tr(:,i) = squeeze(sum(temp_seg(5:15,5:15,:),[1 2]));
end

options.Upper = [max(max(inten_tr)) 600 765 max(max(inten_tr))];
z = -800 : 400 : 800;
z = z';
RMSE_vec = zeros(num_seg,1);
Rsq_vec = zeros(num_seg,1);
mu_vec = zeros(num_seg,1);
On_config = zeros(num_seg,1);
seg_config_RMSE = zeros(num_seg,1);
seg_config_Rsq = zeros(num_seg,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each case represents different segment configuration
% For example case 5 means only the first 3 slices in the segment are fitting
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
% idx_mat{7,1} = logical([0 0 0 1 1]);
% idx_mat{8,1} = logical([0 0 1 1 0]);
% idx_mat{9,1} = logical([0 1 1 0 0]);
% idx_mat{10,1} = logical([1 1 0 0 0]);
z_mat = cell(6,1);
for i = 1:6
    z_mat{i,1} = z(idx_mat{i,1});
end
tic
parfor i = 1 : num_seg  % 16796
    temp = inten_tr(:,i);
    % temp = temp / max(temp);% mustn't normalize
    temp_rmse = nan(6,1);
    temp_Rsq = nan(6,1);
    temp_mu = nan(6,1);
    data_mat = cell(6,1);
    copy_z_mat = z_mat;
    copy_idx_mat = idx_mat;
    % figure % test
    % hold on % test
    % plot(copy_z_mat{1,1},temp) % test 
    z_fine = -800 : 10 : 800;
    for j = 1:6
        data_mat{j,1} = temp(copy_idx_mat{j,1});
        try
            [fit_para_lat1 fg] = fit(copy_z_mat{j,1},data_mat{j,1},fit_model,options);
            z_fit = fit_para_lat1.a1*exp(-((copy_z_mat{j,1}-fit_para_lat1.b1)/fit_para_lat1.c1).^2)+fit_para_lat1.d1;
            temp_rmse(j) = sqrt(sum((z_fit - data_mat{j,1}).^2)/length(z_fit))/(fit_para_lat1.a1+fit_para_lat1.d1);
            temp_Rsq(j) = fg.rsquare;
            temp_mu(j) = fit_para_lat1.b1;
            % z_fit = fit_para_lat1.a1*exp(-((z_fine-fit_para_lat1.b1)/fit_para_lat1.c1).^2)+fit_para_lat1.d1; % test
            % plot(z_fine,z_fit) % test
        catch
        end
    end
    idx_rmse = find(temp_rmse == min(temp_rmse));
    idx_rsq = find(temp_Rsq == max(temp_Rsq));
    mu_vec(i) = temp_mu(idx_rmse(1));
    RMSE_vec(i) = temp_rmse(idx_rmse(1));
    Rsq_vec(i) = temp_Rsq(idx_rsq(1));
    seg_config_RMSE(i) = idx_rmse(1);
    seg_config_Rsq(i) = idx_rsq(1);
end
toc
%%%%%%%%%%%%%%%%%
bin_edge_rmse = linspace(0,0.1,100);
bin_edge_rsq = linspace(0,1,100);
figure
histogram(RMSE_vec,bin_edge_rmse)
xlabel('Normalized RMSE')
figure
histogram(Rsq_vec,bin_edge_rsq)
xlabel('R-squared')
counter_templete = [1 1 1 1 1];
slice_counter_rmse = zeros(num_seg,1);
slice_counter_rsq = zeros(num_seg,1);
for i = 1:num_seg
    slice_counter_rmse(i) = sum(counter_templete(idx_mat{seg_config_RMSE(i),1}));
    slice_counter_rsq(i) = sum(counter_templete(idx_mat{seg_config_Rsq(i),1}));
end
bin_edge_slice = 0.5:1:5.5;

figure
h = histogram(slice_counter_rmse,bin_edge_slice);
% figure
% bar(h.Values)
% set(gca,'YScale','log')
xlabel('Slices per segment')
title('Sorted by rmse')

figure
histogram(slice_counter_rsq,bin_edge_slice);
set(gca, 'YScale', 'log')
xlabel('Slices per segment')
title('Sorted by r-squred')

bin_edge_slice_config = 0.5:1:10.5;
figure
histogram(seg_config_RMSE,bin_edge_slice_config);
save(strcat(data_path,'\seg_config_RMSE_optimal_SM_',scan_mode,'.mat'),'seg_config_RMSE');
% save('seg_config_RMSE_optimal_SM1.mat','seg_config_RMSE')
% save('mu_vec_SM1.mat','mu_vec');