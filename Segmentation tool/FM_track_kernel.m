function PSFdata = FM_track_kernel(img_raw,para,last_seg)
% Gx_oprt=[-1 0 1];
% Gy_oprt=[-1;0;1];
% Gz_oprt=[-1;0;1];
% filter_Gau=para.filter_Gau;
% stepsize=para.stepsize;  % nm
Hbox_size_axi = para.Hbox_size_axi;
Hbox_size_lat = para.Hbox_size_search;
roi_z = para.height_seg;
img_median=zeros(size(img_raw));
width_seg = para.width_seg;
height_seg = para.height_seg;
for i=1:size(img_raw,3)
    img_median(:,:,i)=medfilt2(img_raw(:,:,i));
end
% img_filter=imfilter(img_median,filter_Gau);
% img_filter_div_=double(div_img(img_filter,Gx_oprt,Gy_oprt,Gz_oprt,stepsize));
% img_filter_div_=imfilter(img_filter_div_,filter_Gau);

%% update box center position
last_seg_pos = last_seg.PosSet;
last_seg_vol = last_seg.VolumeSet;
PosSet_update = zeros(size(last_seg_pos));
VolumeSet = cell(size(last_seg_pos,1),1);
NA_idx = [];
for i = 1:size(last_seg_pos,1)
    pos_old_temp = last_seg_pos(i,:);
    try
        box_temp = img_median(pos_old_temp(1)-Hbox_size_lat:pos_old_temp(1)+Hbox_size_lat,pos_old_temp(2)-Hbox_size_lat:pos_old_temp(2)+Hbox_size_lat,pos_old_temp(3)-Hbox_size_axi:pos_old_temp(3)+Hbox_size_axi);
        arr_temp = reshape(box_temp,[((Hbox_size_lat)*2+1)^2*(Hbox_size_axi*2+1) 1]);
        idx_temp = find(arr_temp==max(arr_temp));
        idx_temp = idx_temp(1);
        z_loc = ceil(idx_temp/((Hbox_size_lat)*2+1)^2);
        temp_xy = box_temp(:,:,z_loc);
        [x_loc y_loc] = find(temp_xy==max(max(temp_xy)));
        x_loc = x_loc(1);
        y_loc = y_loc(1);
        x_update = pos_old_temp(1)+(x_loc-(Hbox_size_lat+1));
        y_update = pos_old_temp(2)+(y_loc-(Hbox_size_lat+1));
        z_update = pos_old_temp(3)+(z_loc-(Hbox_size_axi+1));
        PosSet_update(i,:) = [x_update y_update z_update];
        
        VolumeSet{i,1} = img_raw(x_update-width_seg:x_update+width_seg,y_update-width_seg:y_update+width_seg,z_update-height_seg:z_update+height_seg);
        if sum((img_raw(x_update-4:x_update+4,y_update-4:y_update+4,z_update-1:z_update+1)-100),'all')<0.3*sum((last_seg_vol{i,1}(Hbox_size_lat-4:Hbox_size_lat+4,Hbox_size_lat-4:Hbox_size_lat+4,roi_z:roi_z+2)-100),'all')
            NA_idx = [NA_idx i];
        end
    catch
        NA_idx = [NA_idx i];
    end
end
VolumeSet(NA_idx,:) = [];
PosSet_update(NA_idx,:) = [];
FWHMset = zeros(size(PosSet_update));
Manually_sele_idx = [];
IMG_MIP=uint16(img_median-median(img_median,'all'));
xyMIP = max(IMG_MIP,[],3);
PSFdata.VolumeSet = VolumeSet;
PSFdata.xyMIP = xyMIP;
PSFdata.PosSet = PosSet_update;
PSFdata.FWHMset = FWHMset;
PSFdata.Manually_sele_idx = Manually_sele_idx;
end