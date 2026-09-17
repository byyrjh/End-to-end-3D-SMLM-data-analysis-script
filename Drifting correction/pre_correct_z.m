function raw_data = pre_correct_z(raw_data,setting_para)
load(strcat(setting_para.data_path,'segment_data\SM_os_AO_4D_map_m0.mat'));
SlicesPerVol = setting_para.SlicesPerVol;
step_size = setting_para.step_size;
[xyz t_range] = size(LS_os_AO_4D_map);
xy = sqrt(xyz/SlicesPerVol);
LS_os_map_ref = zeros(xy*xy*SlicesPerVol,1);
xybinsize_SM = setting_para.xybinsize;
pixel_size_xy = setting_para.pixel_size_xy;
for i = 1:xy*xy*SlicesPerVol
    counter_t = 0;
    while (counter_t<t_range)
        counter_t = counter_t+1;
        if ~isnan(LS_os_AO_4D_map(i,counter_t))
            LS_os_map_ref(i) = LS_os_AO_4D_map(i,counter_t);
            break;
        end
    end
end


%%%%%%%%%%%  this is 2D case
%{  
for i = 1:t_range
    LS_os_map(:,:,i) = reshape(LS_os_AO_map(:,i),[xy xy]);
end
LS_os_map_ref = zeros(xy,xy);
for x = 1:xy
    for y = 1:xy
        counter_t = 0;
        while (counter_t<t_range)
            counter_t = counter_t+1;
            if ~isnan(LS_os_map(x,y,counter_t))
                LS_os_map_ref(x,y) = LS_os_map(x,y,counter_t);
                break;
            end
        end
    end
end
%}
for i = 1:length(raw_data.z)
    cur_x_bin = ceil(round(raw_data.x(i) / pixel_size_xy) / (xybinsize_SM / pixel_size_xy));
    cur_y_bin = ceil(round(raw_data.y(i) / pixel_size_xy) / (xybinsize_SM / pixel_size_xy));
    cur_z = ceil(raw_data.z(i)/step_size);
    cur_z = min(cur_z,SlicesPerVol);
    cur_t = raw_data.t(i);
    idx_sp = (cur_z-1)*xy*xy + (cur_y_bin-1)*xy + cur_x_bin;
try
    temp_ls_os_m0 = LS_os_AO_4D_map(idx_sp,cur_t)-LS_os_map_ref(idx_sp);
catch
    a = 0;
end
    if ~isnan(temp_ls_os_m0)
        raw_data.z(i) = raw_data.z(i) + temp_ls_os_m0*setting_para.pixel_size_z;
    end
    
%     cur_x_bin = ceil(raw_data.x(i) / pixel_size_xy / (xybinsize_SM / pixel_size_xy));
%     cur_y_bin = ceil(raw_data.y(i) / pixel_size_xy / (xybinsize_SM / pixel_size_xy));
%     cur_t = raw_data.t(i);
%     temp_ls_os_m0 = LS_os_map(cur_x_bin,cur_y_bin,cur_t)-LS_os_map_ref(cur_x_bin,cur_y_bin);
%     if ~isnan(temp_ls_os_m0)
%         raw_data.z(i) = raw_data.z(i) + temp_ls_os_m0*setting_para.pixel_size_z;
%     end

end
end