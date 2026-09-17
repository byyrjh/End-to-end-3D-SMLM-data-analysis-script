function out = fit_interface(fitting_info,fitting_data,config_intm)
tot_len = length(config_intm.map_ptr_x);
fitting_para_tot = config_intm.fitting_para;
fit_slice = config_intm.slice;
fitting_res = zeros(6,tot_len);
for i = 1:ceil(tot_len/4000)
    idx_ini = (i-1)*4000+1;
    idx_end = min(i*4000,tot_len);
    fitting_para = fitting_para_tot((idx_ini-1)*6+1:idx_end*6);
    fitting_para = single(fitting_para);
    fitting_data.fitting_para = fitting_para;
    fitting_data.map_ptr_x = single(config_intm.map_ptr_x(idx_ini:idx_end));
    fitting_data.map_ptr_y = single(config_intm.map_ptr_y(idx_ini:idx_end));
    fitting_data.seg_data = single(config_intm.data(:,:,(idx_ini-1)*fit_slice+1:idx_end*fit_slice));
    fitting_info.slice_num = config_intm.slice;
    fitting_info.fit_offset = config_intm.fit_os;
    fitting_info.seg_num = idx_end-idx_ini+1;% Not more than 4000
    res_temp = GPUfit(fitting_info,fitting_data);
    fitting_res(:,idx_ini:idx_end) = reshape(res_temp.fitting_para,[6 idx_end-idx_ini+1]);
end
out = fitting_res;
end