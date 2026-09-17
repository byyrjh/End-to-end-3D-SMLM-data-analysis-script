function emitter=loc_core_3D(row,col,page,peak_value,aspect_ratio,num_px_thres,ADU_pre_Hmax,Hbox_size,Hbox_size_axi,sim_median_dist)
emitter=[];
global img3D_proc;
while ~isempty(peak_value)
    [cur_peak pos]=max(peak_value);
    row_temp=row(pos);
    col_temp=col(pos);
    page_temp=page(pos);
    row(pos)=[];
    col(pos)=[];
    page(pos)=[];
    peak_value(pos)=[];
    % assume axial resolution is 2.5 times lower than lateral, squeeze PSF
    % in axial to make it isotropic
    vec_dist=sqrt((row-row_temp).^2+(col-col_temp).^2+((page-page_temp)*aspect_ratio/2.5).^2);
    idx=find(vec_dist<Hbox_size);
    vec_dist_psf=vec_dist(idx);
    row_temp=row_temp+Hbox_size;
    col_temp=col_temp+Hbox_size;
    page_temp=page_temp+Hbox_size_axi;
    if (length(idx)>=num_px_thres)&&(median(vec_dist_psf)<sim_median_dist^2)
        img3D_proc(row_temp-Hbox_size:row_temp+Hbox_size,col_temp-Hbox_size:col_temp+Hbox_size,page_temp-Hbox_size_axi:page_temp+Hbox_size_axi)=0;
        emitter=[emitter;row_temp col_temp page_temp double(cur_peak)];
    elseif cur_peak/2==ADU_pre_Hmax
        img3D_proc(row_temp:row_temp,col_temp:col_temp,page_temp:page_temp)=img3D_proc(row_temp:row_temp,col_temp:col_temp,page_temp:page_temp)/3*2;
    end
    row(idx)=[];
    col(idx)=[];
    page(idx)=[];
    peak_value(idx)=[];
end
end