clear;clc;
para.ADU_diff = 0;
para.width_in = 9;
para.height_in = 2;
para.width_search_in = 9;
para.height_search_in = 3;
para.ADU_thres_in = 50;
para.Std_filter_in = 1;
para.Hbox_size_axi = 3;
para.Hbox_size_search = 9;
para.width_seg = 9;
para.height_seg = 2;
para.datatype = 'FM';
path = 'M:\Hao\2023\5_15\2023-05-15_14-47-29_FOV2003\test_m=0_t=0015_det2.tif';
% imgobj = Tiff(path_ary{1,1},'r');
stack_num=100;
slicepS=50;
% for j=slicepS*5+1:slicepS*5+slicepS % kickout first 5 stacks
%     imgobj.setDirectory(j)
%     seed_stack(:,:,j-slicepS*5)=imgobj.read;
% end
% seed_stack = double(seed_stack);
% seed_seg=PSFseg_kernel(seed_stack,para);
load('PosSet.mat');
load('VolumeSet.mat');
last_seg.PosSet = PosSet;
last_seg.VolumeSet = VolumeSet;
imgobj = Tiff(path,'r');
for j=1:stack_num*slicepS
    imgobj.setDirectory(j)
    img_full(:,:,j)=imgobj.read;
end
data_raw = cell(stack_num,1);
data_seg = cell(stack_num,1);
for j = 1:stack_num
    data_raw{j,1} = double(img_full(:,:,(j-1)*slicepS+1:j*slicepS));
end
for j=1:stack_num
    data_seg{j,1} = FM_track_kernel(data_raw{j,1},para,last_seg);
end





