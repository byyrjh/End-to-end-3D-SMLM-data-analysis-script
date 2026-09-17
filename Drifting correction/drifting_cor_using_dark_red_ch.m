clear;clc;close all;
smooth_f=30;% localization precision is independent of smooth binning size
correction_method = 'linear'; % alternative 'mean' 'weighted' 'linear'
outlier_thres_xy=120;
outlier_thres_z=120;
pixel_size_z=10;  % nm
pixel_size_xy=100;  % nm
step_size=400;  % nm
exposuretime=30;  % ms
num_photon_factor=30;
slice_per_stack=25;
stationary_pos = 5;
vol_per_hyper = 200;
FM_outlier = [];
origin_drifting=1;
scan_dir = 'm0';
x_shift = -550;
y_shift = 170;
search_margin = 300; % nm
%%
datasave_path='M:\Hao\2022\2022_12_11\cell1\segment_data\';
load(strcat(datasave_path,'map_ptr_t_FM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_x_FM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_y_FM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_z_FM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_t_SM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_x_SM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_y_SM_',scan_dir,'.mat'));
load(strcat(datasave_path,'map_ptr_z_SM_',scan_dir,'.mat'));
load(strcat(datasave_path,'fitting_result_FM_',scan_dir,'.mat'));
%% fiducial marker track
x_pos=(map_ptr_x_FM+fitting_results(1,:)')*pixel_size_xy;
y_pos=(map_ptr_y_FM+fitting_results(2,:)')*pixel_size_xy;
z_pos=map_ptr_z_FM*step_size+fitting_results(3,:)'*pixel_size_z;
cor_lab=[x_pos y_pos z_pos];
data_full=[cor_lab fitting_results(4,:)'*num_photon_factor fitting_results(6,:)'*pixel_size_z map_ptr_t_FM];
data_full=double(round(data_full));
maxdisp=500;  % nm
param=struct('mem',0,'dim',3,'good',1,'quiet',1);
track_output_m0=track(data_full,maxdisp,param);
traj_m0=zeros(max(track_output_m0(:,7)),1);
for i=1:size(track_output_m0,1)
    traj_m0(track_output_m0(i,7))=traj_m0(track_output_m0(i,7))+1;
end
marker_idx_m0=find(traj_m0==max(map_ptr_t_FM));
origin_arr=zeros(length(marker_idx_m0),3);
%% raw trajactory generation
trace_x_m0=[];
trace_y_m0=[];
trace_z_m0=[];
trace_inten_m0=[];
for i=1:length(marker_idx_m0)
    temp_x=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),1);
    temp_y=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),2);
    temp_z=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),3);
    temp_inten=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),4);
    trace_x_m0=[trace_x_m0 temp_x];
    trace_y_m0=[trace_y_m0 temp_y];
    trace_z_m0=[trace_z_m0 temp_z];
    origin_arr(i,1)=mean(temp_x(6:30));
    origin_arr(i,2)=mean(temp_y(6:30));
    origin_arr(i,3)=mean(temp_z(6:30));
    trace_inten_m0=[trace_inten_m0 temp_inten];
end
num_vol=length(trace_inten_m0);
vol_exposure_t=exposuretime*slice_per_stack*2/1000;
duration_meas=vol_exposure_t*num_vol/60; % min
t_axis=linspace(0,duration_meas-exposuretime*slice_per_stack*2/1000/60,num_vol);
figure
hold on
xlabel('Time / min')
ylabel('Drifting / nm')
box on 
grid on
title('Drifting in x raw')
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_x_m0(:,i))
end

figure
hold on
xlabel('Time / min')
ylabel('Drifting / nm')
box on 
grid on
title('Drifting in y raw')
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_y_m0(:,i))
end

figure
hold on
xlabel('Time / min')
ylabel('Drifting / nm')
box on 
grid on
title('Drifting in z raw')
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_z_m0(:,i))
end

figure
hold on
title('Photobleaching')
xlabel('Time / min')
ylabel('Number of photons')
box on 
grid on
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,trace_inten_m0(:,i))
end
%% light sheet drifting trajactory
LS_drifting=zeros(max(map_ptr_t_FM),length(marker_idx_m0));
load(strcat(datasave_path,'LS_os_',scan_dir,'.mat'));
for i=1:length(marker_idx_m0)
    LS_drifting(:,i)=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),5);
end
figure
hold on
title('Light sheet offset')
xlabel('Time / min')
ylabel('Distance / nm')
box on 
grid on
set(gca,'FontSize',20)
LS_drifting = data_full(:,5);
LS_drifting = reshape(LS_drifting,[length(LS_drifting)/num_vol num_vol])';
for i=1:length(marker_idx_m0)
    plot(t_axis,LS_os(i,:)*10);
%     plot(t_axis,LS_drifting(:,i));
end

LS_os_matlab=zeros(size(map_ptr_t_SM));
for i=1:length(map_ptr_t_SM)
    vol_idx_LS=map_ptr_t_SM(i);
    x=map_ptr_x_SM(i);
    y=map_ptr_y_SM(i);
    x_FM=map_ptr_x_FM((vol_idx_LS-1)*length(marker_idx_m0)+1:vol_idx_LS*length(marker_idx_m0));
    y_FM=map_ptr_y_FM((vol_idx_LS-1)*length(marker_idx_m0)+1:vol_idx_LS*length(marker_idx_m0));
    dist=sqrt((x-x_FM).^2+(y-y_FM).^2);
    cur_LS_os=LS_drifting(vol_idx_LS,:);
    LS_os_matlab(i)=sum(1./dist.*cur_LS_os'/sum(1./dist));
end

%% smoothed trajectory generation 
vol_per_hyper = vol_per_hyper-stationary_pos;
num_hyperstack = num_vol/vol_per_hyper;
drifting_corr_x = zeros(num_vol,length(marker_idx_m0));
drifting_corr_y = zeros(num_vol,length(marker_idx_m0));
drifting_corr_z = zeros(num_vol,length(marker_idx_m0));
for i=1:length(marker_idx_m0)
    for j = 1:num_hyperstack
        drifting_corr_x(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i)=smooth(trace_x_m0(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i),smooth_f);
        drifting_corr_y(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i)=smooth(trace_y_m0(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i),smooth_f);
        drifting_corr_z(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i)=smooth(trace_z_m0(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i),smooth_f);
    end
end

figure
hold on
grid on
box on
set(gca,'FontSize',20)
title('Drifting in x (mean)')
for i=1:length(marker_idx_m0)
    plot(t_axis,drifting_corr_x(:,i))
end
xlabel('Time / min')
ylabel('Drifting / nm')

figure
hold on
grid on
box on
set(gca,'FontSize',20)
title('Drifting in y (mean)')
for i=1:length(marker_idx_m0)
    plot(t_axis,drifting_corr_y(:,i))
end
xlabel('Time / min')
ylabel('Drifting / nm')

figure
hold on
grid on
box on
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,drifting_corr_z(:,i))
end
xlabel('Time / min')
ylabel('Drifting / nm')
title('Drifting in z (mean)')

figure
hold on
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
scatter3(origin_arr(i,1)/1000,origin_arr(i,2)/1000,origin_arr(i,3)/1000,'LineWidth',4);
end
axis equal
xlabel('x / µm')
ylabel('y / µm')
zlabel('z / µm')
box on
grid on
title('Fiducial marker distribution')

%% localization precision generation
trace_x_corrected=trace_x_m0-drifting_corr_x; % trace_x_m0 is raw trajectory
trace_y_corrected=trace_y_m0-drifting_corr_y;
trace_z_corrected=trace_z_m0-drifting_corr_z;

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

for i=1:length(marker_idx_m0)
    for j = 1:num_hyperstack
        drifting_corr_x(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i)=smooth(trace_x_m0(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i),smooth_f);
        drifting_corr_y(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i)=smooth(trace_y_m0(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i),smooth_f);
        drifting_corr_z(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i)=smooth(trace_z_m0(1+(j-1)*vol_per_hyper:j*vol_per_hyper,i),smooth_f);
    end
end

trace_x_corrected=trace_x_m0-drifting_corr_x+ones(num_vol,1)*drifting_corr_x(1,:);
trace_y_corrected=trace_y_m0-drifting_corr_y+ones(num_vol,1)*drifting_corr_y(1,:);
trace_z_corrected=trace_z_m0-drifting_corr_z+ones(num_vol,1)*drifting_corr_z(1,:);
figure
hold on
xlabel('Time / min')
ylabel('Drifting / nm')
box on 
grid on
title('Drifting data for correction')
set(gca,'FontSize',20)
for i=1:length(marker_idx_m0)
    plot(t_axis,drifting_corr_z(:,i)-drifting_corr_z(1,i))
end
%%
%%%%%%%%%%%%%%%%%%%   correct for single molecules   %%%%%%%%%%%%%%%%%%%%%%
load(strcat(datasave_path,'fitting_result_SM_',scan_dir,'.mat'));
load(strcat(datasave_path,'ChiSq_SM_',scan_dir,'.mat'));
x_pos_SM=double((map_ptr_x_SM+fitting_results(1,:)')*pixel_size_xy);
y_pos_SM=double((map_ptr_y_SM+fitting_results(2,:)')*pixel_size_xy);
z_pos_SM=double(map_ptr_z_SM*step_size+fitting_results(3,:)'*pixel_size_z);
inten_SM=double(fitting_results(4,:)'*num_photon_factor);
cor_lab_SM=[x_pos_SM y_pos_SM z_pos_SM];
data_full=[cor_lab_SM inten_SM fitting_results(6,:)'*pixel_size_z map_ptr_t_SM];
for i=1:length(marker_idx_m0)
    search_cent_x = drifting_corr_x(1,i)+x_shift;
    search_cent_y = drifting_corr_y(1,i)+y_shift;
    search_cent_z = drifting_corr_z(1,i);
    search_range_x = abs(drifting_corr_x(num_vol,i)-drifting_corr_x(1,i))+search_margin;
    search_range_y = abs(drifting_corr_y(num_vol,i)-drifting_corr_y(1,i))+search_margin;
    search_range_z = abs(drifting_corr_z(num_vol,i)-drifting_corr_z(1,i))+search_margin;
    idx_x = (x_pos_SM>search_cent_x-search_range_x)&(x_pos_SM<search_cent_x+search_range_x);
    idx_y = (y_pos_SM>search_cent_y-search_range_y)&(y_pos_SM<search_cent_y+search_range_y);
    idx_z = (z_pos_SM>search_cent_z-search_range_z)&(z_pos_SM<search_cent_z+search_range_z);
    idx_tot = idx_x|idx_y|idx_z;
    cur_sub_vol = double(round(data_full(idx_tot,:)));
    track_output_m0=track(cur_sub_vol,maxdisp,param);
    traj_m0=zeros(max(track_output_m0(:,7)),1);
    for j=1:size(track_output_m0,1)
        traj_m0(track_output_m0(j,7))=traj_m0(track_output_m0(j,7))+1;
    end
    marker_idx_m0=find(traj_m0==max(map_ptr_t_FM));

end



origin_arr=zeros(length(marker_idx_m0),3);




photon_bin_size=0:50:4000;
figure 
hold on
histogram(inten_SM,photon_bin_size)
xlim([0 4000])
xlabel('Number of photons')
ylabel('Occurrence')
box on
set(gca,'FontSize',20)
t_SM=double(map_ptr_t_SM);
nan_x=isnan(x_pos_SM);
nan_y=isnan(y_pos_SM);
nan_z=isnan(z_pos_SM);
x_outlier_idx = (fitting_results(1,:)>4)';% single molecule is found at the edge of segment
y_outlier_idx = (fitting_results(2,:)>4)';% single molecule is found at the edge of segment
z_outlier_idx = (fitting_results(3,:)>40)';% single molecule is found far away from center slice of segment
os_outlier_idx = (fitting_results(6,:)>100)';
kick_out_idx=nan_x|nan_y|nan_z|x_outlier_idx|y_outlier_idx|z_outlier_idx|os_outlier_idx;
x_pos_SM(kick_out_idx)=[];
y_pos_SM(kick_out_idx)=[];
z_pos_SM(kick_out_idx)=[];
inten_SM(kick_out_idx)=[];
inten_SM=ones(size(inten_SM));
t_SM(kick_out_idx)=[];
drifting_corr_x(:,FM_outlier) = [];
drifting_corr_y(:,FM_outlier) = [];
drifting_corr_z(:,FM_outlier) = [];
%% correction kernel
switch correction_method
    case 'weighted'
        for i=1:num_vol
            idx=(t_SM==i);
            SM_x=x_pos_SM(idx);
            SM_y=y_pos_SM(idx);
            SM_z=z_pos_SM(idx);
            FM_x_abs=drifting_corr_x(i,:);
            FM_y_abs=drifting_corr_y(i,:);
            FM_z_abs=drifting_corr_z(i,:);
            FM_x_rel=drifting_corr_x(i,:)-drifting_corr_x(origin_drifting,:);
            FM_y_rel=drifting_corr_y(i,:)-drifting_corr_y(origin_drifting,:);
            FM_z_rel=drifting_corr_z(i,:)-drifting_corr_z(origin_drifting,:);
            dist_x=SM_x-ones(size(SM_x))*FM_x_abs;
            dist_y=SM_y-ones(size(SM_y))*FM_y_abs;
            dist_z=SM_z-ones(size(SM_z))*FM_z_abs;
            mat_dist=sqrt(dist_x.^2+dist_y.^2+dist_z.^2);
            mat_weight=(1./mat_dist)./(sum(1./mat_dist,2)*ones(size(FM_x_abs)));
            corr_x=sum(ones(size(SM_x))*FM_x_rel.*mat_weight,2);
            corr_y=sum(ones(size(SM_y))*FM_y_rel.*mat_weight,2);
            corr_z=sum(ones(size(SM_z))*FM_z_rel.*mat_weight,2);
            x_pos_SM(idx)=x_pos_SM(idx)-corr_x;
            y_pos_SM(idx)=y_pos_SM(idx)-corr_y;
            z_pos_SM(idx)=z_pos_SM(idx)-corr_z;
        end
    case 'mean'
        for i=1:num_vol
            idx=(t_SM==i);
            SM_x=x_pos_SM(idx);
            SM_y=y_pos_SM(idx);
            SM_z=z_pos_SM(idx);
            FM_x_rel=drifting_corr_x(i,:)-drifting_corr_x(origin_drifting,:);
            FM_y_rel=drifting_corr_y(i,:)-drifting_corr_y(origin_drifting,:);
            FM_z_rel=drifting_corr_z(i,:)-drifting_corr_z(origin_drifting,:);
            x_pos_SM(idx)=x_pos_SM(idx)-mean(FM_x_rel);
            y_pos_SM(idx)=y_pos_SM(idx)-mean(FM_y_rel);
            z_pos_SM(idx)=z_pos_SM(idx)-mean(FM_z_rel);
        end
    case 'linear'
        for i=1:num_vol
            idx=(t_SM==i); 
            axis_z_cur=drifting_corr_z(i,:)';
            [axis_z_cur,rank_order]=sort(axis_z_cur);
            drift_x_cur=drifting_corr_x(i,:)'-drifting_corr_x(origin_drifting,:)';
            drift_y_cur=drifting_corr_y(i,:)'-drifting_corr_y(origin_drifting,:)';
            drift_x_cur=drift_x_cur(rank_order);
            drift_y_cur=drift_y_cur(rank_order);
            coef_x = polyfit(axis_z_cur,drift_x_cur,1);
            coef_y = polyfit(axis_z_cur,drift_y_cur,1);
            FM_z_rel=drifting_corr_z(i,:)'-drifting_corr_z(origin_drifting,:)';
            FM_z_rel = FM_z_rel(rank_order);
            x_pos_SM(idx)=x_pos_SM(idx)-(z_pos_SM(idx)*coef_x(1)+coef_x(2));
            y_pos_SM(idx)=y_pos_SM(idx)-(z_pos_SM(idx)*coef_y(1)+coef_y(2));
            z_pos_SM(idx)=z_pos_SM(idx)-mean(FM_z_rel);           
        end

end


FM_x=reshape(trace_x_corrected,[num_vol*length(marker_idx_m0) 1]);
FM_y=reshape(trace_y_corrected,[num_vol*length(marker_idx_m0) 1]);
FM_z=reshape(trace_z_corrected,[num_vol*length(marker_idx_m0) 1]);
FM_inten=10*ones(size(FM_z));
FM_t=1*ones(size(FM_z));

theta = 33/180*pi;
t = [ 1          0            0;
    0      cos(theta)  -sin(theta);
    0      sin(theta)   cos(theta) ];
full_SM_data=[x_pos_SM y_pos_SM z_pos_SM inten_SM t_SM];
full_FM_data=[FM_x FM_y FM_z FM_inten FM_t];
% full_SM_data(:,1:3)=(t*full_SM_data(:,1:3)')';
% full_FM_data(:,1:3)=(t*full_FM_data(:,1:3)')';
full_SM_data=round(full_SM_data);
full_FM_data=round(full_FM_data);
% xmin = 118;
% xmax = 21768;
% ymin = -6307;
% ymax = 16254;
% zmin = 9491;
% zmax = 17669;
% ROI_idx = (full_SM_data(:,1)>xmin) & (full_SM_data(:,1)<xmax) & full_SM_data(:,2)>ymin & full_SM_data(:,2)<ymax & full_SM_data(:,3)>zmin & full_SM_data(:,3)<zmax ;
% full_SM_data = full_SM_data(ROI_idx,:);
%%
full_SM_data=[full_SM_data; full_FM_data];
tg=fopen(strcat(datasave_path,'SM_loc_',scan_dir,correction_method,'_smoothF_lab_cor_',num2str(smooth_f),'.txt'),'wt');
for j=1:size(full_SM_data,1)
    fprintf(tg,'%f\t%f\t%f\t%f\t%f\n',full_SM_data(j,1),full_SM_data(j,2),full_SM_data(j,3),full_SM_data(j,4),full_SM_data(j,5));
end
fclose(tg);
