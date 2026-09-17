function [FM_trace]= get_FM_trace(dir, scan, setting_para)
%%
% obtain oringinal fiducial marker trajectory, number of photons trance and
% light sheet offset 
% outlier spike is replaced with average position of its surroundings
%%
smooth_f = setting_para.smooth_f;% localization precision is independent of smooth binning size
outlier_thres_xy = setting_para.outlier_thres_xy;
outlier_thres_z = setting_para.outlier_thres_z;
pixel_size_z = setting_para.pixel_size_z;  % nm
pixel_size_xy = setting_para.pixel_size_xy;  % nm
step_size = setting_para.step_size;  % nm
num_photon_factor = setting_para.num_photon_factor;
switch scan
    case 'm0'
        FM_outlier = setting_para.FM_outlier{1,1};
    case 'm1'
        FM_outlier = setting_para.FM_outlier{2,1};
end

load(strcat(dir,'map_ptr_t_FM_',scan,'.mat'));
load(strcat(dir,'map_ptr_x_FM_',scan,'.mat'));
load(strcat(dir,'map_ptr_y_FM_',scan,'.mat'));
load(strcat(dir,'map_ptr_z_FM_',scan,'.mat'));
load(strcat(dir,'fitting_result_FM_',scan,'.mat'));
load(strcat(dir,'fitting_info_',scan,'.mat'));
x_pos=double(round((map_ptr_x_FM+fitting_results(1,:)')*pixel_size_xy));
y_pos=double(round((map_ptr_y_FM+fitting_results(2,:)')*pixel_size_xy));
z_pos=double(round(map_ptr_z_FM*step_size+fitting_results(3,:)'*pixel_size_z));
inten = double(round(fitting_results(4,:)'*num_photon_factor));
LS_os = double(round(fitting_results(6,:)'*pixel_size_z)); % smoothed in c code
% num_FM_tr = num_FM/num_vol;
% x_pos = reshape(x_pos,[num_FM_tr num_vol])';
% y_pos = reshape(y_pos,[num_FM_tr num_vol])';
% z_pos = reshape(z_pos,[num_FM_tr num_vol])';

%% 
cor_lab=[x_pos y_pos z_pos];
data_full=[cor_lab inten LS_os map_ptr_t_FM];
data_full=double(round(data_full));
maxdisp=1000;  % nm
param=struct('mem',0,'dim',3,'good',1,'quiet',1);
track_output_m0=track(data_full,maxdisp,param);
traj_m0=zeros(max(track_output_m0(:,7)),1);
for i=1:size(track_output_m0,1)
    traj_m0(track_output_m0(i,7))=traj_m0(track_output_m0(i,7))+1;
end
marker_idx_m0=find(traj_m0==max(map_ptr_t_FM));
trace_x_m0=[];
trace_y_m0=[];
trace_z_m0=[];
trace_inten_m0=[];
trace_os_m0 = [];
for i=1:length(marker_idx_m0)
    temp_x=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),1);
    temp_y=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),2);
    temp_z=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),3);
    temp_inten=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),4);
    temp_os=track_output_m0(track_output_m0(:,7)==marker_idx_m0(i),5);
    trace_x_m0=[trace_x_m0 temp_x];
    trace_y_m0=[trace_y_m0 temp_y];
    trace_z_m0=[trace_z_m0 temp_z];
    trace_inten_m0=[trace_inten_m0 temp_inten];
    trace_os_m0 = [trace_os_m0 temp_os];
end
x_pos = trace_x_m0;
y_pos = trace_y_m0;
z_pos = trace_z_m0;
[len_FM_tr, num_FM_tr] = size(x_pos);
%%
vol_per_hyper = setting_para.vol_per_hyper-setting_para.stationary_pos;
drifting_corr_x = smooth_gen(x_pos, vol_per_hyper, setting_para.piecepos, smooth_f);
drifting_corr_y = smooth_gen(y_pos, vol_per_hyper, setting_para.piecepos, smooth_f);
drifting_corr_z = smooth_gen(z_pos, vol_per_hyper, setting_para.piecepos, smooth_f);

trace_x_corrected=x_pos-drifting_corr_x; % x_pos is raw trajectory
trace_y_corrected=y_pos-drifting_corr_y;
trace_z_corrected=z_pos-drifting_corr_z;
for i=1:num_FM_tr   % kick out localization spikes
    idx_temp=(abs(trace_x_corrected(:,i))>outlier_thres_xy)|(abs(trace_y_corrected(:,i))>outlier_thres_xy)|(abs(trace_z_corrected(:,i))>outlier_thres_z);
    idx_pos=find(idx_temp==true);
    if ~isempty(idx_pos)
        for j=1:length(idx_pos)
            if (idx_pos(j)>smooth_f/2 && (idx_pos(j)+smooth_f/2)<num_vol)
                temp_x=x_pos(idx_pos(j)-smooth_f/2:idx_pos(j)+smooth_f/2,i);
                temp_y=y_pos(idx_pos(j)-smooth_f/2:idx_pos(j)+smooth_f/2,i);
                temp_z=z_pos(idx_pos(j)-smooth_f/2:idx_pos(j)+smooth_f/2,i);
                temp_x(smooth_f/2+1)=[];
                temp_y(smooth_f/2+1)=[];
                temp_z(smooth_f/2+1)=[];
                x_pos(idx_pos(j),i)=mean(temp_x);
                y_pos(idx_pos(j),i)=mean(temp_y);
                z_pos(idx_pos(j),i)=mean(temp_z);
            end
        end
    end
end

x_pos(:,FM_outlier) = [];
y_pos(:,FM_outlier) = [];
z_pos(:,FM_outlier) = [];
trace_inten_m0(:,FM_outlier) = [];
trace_os_m0(:,FM_outlier) = [];

FM_trace.x = x_pos;
FM_trace.y = y_pos;
FM_trace.z = z_pos;
FM_trace.inten = trace_inten_m0;
FM_trace.LS_os = trace_os_m0;
% assumption 1. RI change induced fluctuation of detection focal plane
% cannot be ignored thus has to be corrected using light sheet offset
% fitting
% the origin of light sheet offset is the center of PSF_det
% PSF_exc(z-z^hat-z_c) the light sheet offset z^hat share the same sign of z_c
% z^hat >0 means light sheet offset is in the direction z_c increases
% Increasing z^hat (Delta z^hat = z^hat(t1)-z^hat(t0)>0) implies decreasing z_c due to RI change
% Therefore z_c_corrected = z_c + Delta z^hat

% conundrum: Decomposition of apparent light sheet offset into fluctuation
% and plane components does not mean the mechanical misalignment only contributes
% to plane component though light sheet misalignment can be represented by
% a plane because evolution of RI induced focal plane change very like to
% cause piston and tip-tile aberration that have the same effect as
% mechanical misalignment
end