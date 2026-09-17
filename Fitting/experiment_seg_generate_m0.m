clear;clc;%close all
% output files seg_data(total) map_ptr_x(total) map_ptr_y(total)
% map_ptr_z(total) for m=1 flip seg data flip z integer
% log data col_1 time col_2 m=0 or 1 
% PSFdata is created from 2 hyperstacks m=0 and m=1
% emitter is arranged as m=0 t=0 m=0 t=1...m=0 t=end, m=1 t=0 m=1 t=1...m=1 t=end
%%
scriptFolder = fileparts(mfilename('fullpath'));
packageFolder = fileparts(scriptFolder);
prt_file_path = fullfile(packageFolder,'data\');
data_path=prt_file_path;
prt_file_fullpath=strcat(prt_file_path,'test_prt.txt');
pixelsize=100;  % nm
stepsize=400;  % nm
stationary_pos=1;
%%
fileID = fopen(prt_file_fullpath,'r');
formatSpec = '%s';
prt_content=convertCharsToStrings(fscanf(fileID,formatSpec));
ROI_identifier_start='[L,T,R,B]:';
ROI_identifier_end='Movement';
pt_idf_start=strfind(prt_content,convertCharsToStrings(ROI_identifier_start));
pt_idf_end=strfind(prt_content,convertCharsToStrings(ROI_identifier_end));
ROIinfo=prt_content{1}(pt_idf_start(1)+10:pt_idf_end(1)-1);
seperator=strfind(ROIinfo,',');
L=str2num(ROIinfo(1:seperator(1)-1));
T=str2num(ROIinfo(seperator(1)+1:seperator(2)-1));
R=str2num(ROIinfo(seperator(2)+1:seperator(3)-1));
B=str2num(ROIinfo(seperator(3)+1:length(ROIinfo)));
if ~isfolder(strcat(prt_file_path,'setup_calibration'))
    mkdir(strcat(prt_file_path,'setup_calibration'))
end
if ~isfolder(strcat(prt_file_path,'segment_data'))
    mkdir(strcat(prt_file_path,'segment_data'))
end
load(strcat(prt_file_path,'\gain_det1.mat'));
load(strcat(prt_file_path,'\darkimage_cam1.mat'));
gain=single(sCMOSgain(T:B,L:R));
offset=single(imgseq{1,1}(T:B,L:R));
var=single(imgseq{2,1}(T:B,L:R));
save(strcat(prt_file_path,'setup_calibration\camera1_cali.mat'),'gain','offset','var');
load(strcat(prt_file_path,'\gain_det2.mat'));
load(strcat(prt_file_path,'\darkimage_cam2.mat'));
gain=single(sCMOSgain(T:B,L:R));
offset=single(imgseq{1,1}(T:B,L:R));
var=single(imgseq{2,1}(T:B,L:R));
[cam_size_x cam_size_y] = size(offset);
save(strcat(prt_file_path,'setup_calibration\camera2_cali.mat'),'gain','offset','var');
%% fiducial marker prep
load(strcat(data_path,'PSFdata_FM_m0.mat'));
vol_num=size(data.data,1);
vol_idx=linspace(1,vol_num,vol_num);
vol_per_hyper = size(data.data,1)/size(data.path,1);
spike_idx=(mod(vol_idx,vol_per_hyper)<stationary_pos&mod(vol_idx,vol_per_hyper)>0)';
data.data(spike_idx) = [];
overall_stack_num=size(data.data,1);
SlicesPerSeg=size(data.data{1,1}.VolumeSet{1,1},3);
slice_num_FM = single(SlicesPerSeg);
num_emitter=zeros(overall_stack_num,1);
seg_size=size(data.data{1,1}.VolumeSet{1,1},1);
for stack_idx=1:overall_stack_num
    num_emitter(stack_idx)=size(data.data{stack_idx,1}.VolumeSet,1);
end
num_emitter_tot=sum(num_emitter,'all');
seg_data_FM=zeros(seg_size,seg_size,num_emitter_tot*SlicesPerSeg);
map_ptr_x_FM=zeros(num_emitter_tot,1);
map_ptr_y_FM=zeros(num_emitter_tot,1);
map_ptr_z_FM=zeros(num_emitter_tot,1);
map_idx=linspace(1,num_emitter_tot,num_emitter_tot)';
map_time=zeros(num_emitter_tot,1);
for stack_idx=1:overall_stack_num
    map_ptr_x_FM(sum(num_emitter(1:stack_idx-1,1))+1:sum(num_emitter(1:stack_idx,1)))=data.data{stack_idx,1}.PosSet(:,1);
    map_ptr_y_FM(sum(num_emitter(1:stack_idx-1,1))+1:sum(num_emitter(1:stack_idx,1)))=data.data{stack_idx,1}.PosSet(:,2);
    map_ptr_z_FM(sum(num_emitter(1:stack_idx-1,1))+1:sum(num_emitter(1:stack_idx,1)))=data.data{stack_idx,1}.PosSet(:,3);
    map_time(sum(num_emitter(1:stack_idx-1,1))+1:sum(num_emitter(1:stack_idx,1)))=stack_idx;
    for i=1:num_emitter(stack_idx,1)
        seg_data_FM(:,:,sum(num_emitter(1:stack_idx-1,1))*SlicesPerSeg+(i-1)*SlicesPerSeg+1:sum(num_emitter(1:stack_idx-1,1))*SlicesPerSeg+i*SlicesPerSeg)=data.data{stack_idx,1}.VolumeSet{i,1};
    end
end
%% FM trace
data_full=[map_ptr_x_FM map_ptr_y_FM map_ptr_z_FM*stepsize/pixelsize map_idx map_time];
maxdisp=10;  % lateral pixel
param=struct('mem',0,'dim',3,'good',1,'quiet',1);
track_output=track(data_full,maxdisp,param);
traj=zeros(max(track_output(:,6)),1);
for i=1:size(track_output,1)
    traj(track_output(i,6))=traj(track_output(i,6))+1;
end
marker_idx=find(traj==overall_stack_num);
valid_FM_idx=zeros(num_emitter_tot,1);
for i=1:length(marker_idx)
    valid_FM_idx=valid_FM_idx|(track_output(:,6)==marker_idx(i));
end

valid_FM=track_output(valid_FM_idx,:); % x y z old_idx t trace_idx
[~, b]=sort(valid_FM(:,5));
valid_FM=valid_FM(b,:);
num_valid_FM=length(b);
num_FM=length(marker_idx);
num_t=num_valid_FM/num_FM;
for i=1:num_t
    temp=valid_FM((i-1)*num_FM+1:i*num_FM,:);
    [~, b_idx]=sort(temp(:,6));
    valid_FM((i-1)*num_FM+1:i*num_FM,:)=temp(b_idx,:);
end
% for i=1:num_t
%     temp=valid_FM((i-1)*num_FM+1:i*num_FM,:);
%     [~, b_idx]=sort(temp(:,2));
%     valid_FM((i-1)*num_FM+1:i*num_FM,:)=temp(b_idx,:);
% end

seg_data=zeros(seg_size,seg_size,num_valid_FM*SlicesPerSeg);
for i=1:num_valid_FM
    cur_idx=valid_FM(i,4);
    seg_data(:,:,(i-1)*SlicesPerSeg+1:i*SlicesPerSeg)=seg_data_FM(:,:,(cur_idx-1)*SlicesPerSeg+1:cur_idx*SlicesPerSeg);
end
map_ptr_x_FM=single(valid_FM(:,1)); % FM data array structure 
map_ptr_y_FM=single(valid_FM(:,2)); %[FM_1_t0 FM_2_t0 FM_3_t0...FM_n_t0 FM_1_t1 FM_2_t1 FM_3_t1...FM_n_t1 FM_1_t2 FM_2_t2 FM_3_t2... ]'
map_ptr_z_FM=single(valid_FM(:,3)/stepsize*pixelsize);
map_ptr_t_FM=single(valid_FM(:,5));
seg_data_FM=single(seg_data);
save(strcat(data_path,'segment_data\seg_data_FM_m0.mat'),'seg_data_FM');
save(strcat(data_path,'segment_data\map_ptr_x_FM_m0.mat'),'map_ptr_x_FM');
save(strcat(data_path,'segment_data\map_ptr_y_FM_m0.mat'),'map_ptr_y_FM');
save(strcat(data_path,'segment_data\map_ptr_t_FM_m0.mat'),'map_ptr_t_FM');
save(strcat(data_path,'segment_data\map_ptr_z_FM_m0.mat'),'map_ptr_z_FM');
fprintf('The number of fiducial markers per volume is %d\n',length(marker_idx));
fprintf('The total number of volume is %d\n',overall_stack_num);
fprintf('The total number of FM segments is %d\n',length(map_ptr_x_FM));
num_FM = single(length(map_ptr_x_FM));
num_vol = single(overall_stack_num);
%% single molecule prep
load(strcat(data_path,'PSFdata_SM_m0.mat'));
data.data(spike_idx) = [];
overall_stack_num=size(data.data,1);
SlicesPerSeg=size(data.data{1,1}.VolumeSet{1,1},3);
slice_num_SM = single(SlicesPerSeg);
num_emitter=zeros(overall_stack_num,1);
seg_size=size(data.data{1,1}.VolumeSet{1,1},1);
for stack_idx=1:overall_stack_num
    num_emitter(stack_idx)=size(data.data{stack_idx,1}.VolumeSet,1);
end
num_emitter_tot=sum(num_emitter,'all');
seg_data_SM=zeros(seg_size,seg_size,num_emitter_tot*SlicesPerSeg);
map_ptr_x_SM=zeros(num_emitter_tot,1);
map_ptr_y_SM=zeros(num_emitter_tot,1);
map_ptr_z_SM=zeros(num_emitter_tot,1);
map_ptr_t_SM=zeros(num_emitter_tot,1);
for stack_idx=1:overall_stack_num
    map_ptr_x_SM(sum(num_emitter(1:stack_idx-1,1))+1:sum(num_emitter(1:stack_idx,1)))=data.data{stack_idx,1}.PosSet(:,1);
    map_ptr_y_SM(sum(num_emitter(1:stack_idx-1,1))+1:sum(num_emitter(1:stack_idx,1)))=data.data{stack_idx,1}.PosSet(:,2);
    map_ptr_z_SM(sum(num_emitter(1:stack_idx-1,1))+1:sum(num_emitter(1:stack_idx,1)))=data.data{stack_idx,1}.PosSet(:,3);
    map_ptr_t_SM(sum(num_emitter(1:stack_idx-1,1))+1:sum(num_emitter(1:stack_idx,1)))=stack_idx;
    for i=1:num_emitter(stack_idx,1)
        seg_data_SM(:,:,sum(num_emitter(1:stack_idx-1,1))*SlicesPerSeg+(i-1)*SlicesPerSeg+1:sum(num_emitter(1:stack_idx-1,1))*SlicesPerSeg+i*SlicesPerSeg)=data.data{stack_idx,1}.VolumeSet{i,1};
    end
end
map_ptr_x_SM=single(map_ptr_x_SM);
map_ptr_y_SM=single(map_ptr_y_SM);
map_ptr_z_SM=single(map_ptr_z_SM);
map_ptr_t_SM=single(map_ptr_t_SM);
seg_data_SM=single(seg_data_SM);

save(strcat(data_path,'segment_data\seg_data_SM_m0.mat'),'seg_data_SM', '-v7.3');
save(strcat(data_path,'segment_data\map_ptr_x_SM_m0.mat'),'map_ptr_x_SM');
save(strcat(data_path,'segment_data\map_ptr_y_SM_m0.mat'),'map_ptr_y_SM');
save(strcat(data_path,'segment_data\map_ptr_z_SM_m0.mat'),'map_ptr_z_SM');
save(strcat(data_path,'segment_data\map_ptr_t_SM_m0.mat'),'map_ptr_t_SM');
fprintf('The total number of SM segments is %d\n',length(map_ptr_x_SM));
fprintf('The recommended number of fitting blocks for single molecules is %d \n',ceil(length(map_ptr_x_SM)/5000))
fprintf('The camera size is %d * %d \n',cam_size_x, cam_size_y)

num_SM = single(length(map_ptr_x_SM));
camsize_x = single(cam_size_x);
camsize_y = single(cam_size_y);
stationary_pos = single(stationary_pos);
vol_per_hyper = single(vol_per_hyper);
save(strcat(data_path,'segment_data\fitting_info_m0.mat'),'num_FM','num_SM','slice_num_FM','slice_num_SM','num_vol','camsize_x','camsize_y','stationary_pos','vol_per_hyper');