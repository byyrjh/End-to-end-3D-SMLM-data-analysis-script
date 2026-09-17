function registered_data = register_m1(raw_data, raw_drifting_tr, vol_per_hyper, piecepos, coef_struc)
t_SM = raw_data.t;
num_vol = max(raw_data.t);
drift_z_m0 = raw_drifting_tr.z+raw_drifting_tr.LS_os_rel; % drift in z is corrected in terms of focal plane fluctuation
for i=1:num_vol
    idx=(t_SM==i);
    pos_idx = find(t_SM==i);
    hyper_idx = ceil(i/vol_per_hyper);
    seg_idx = (hyper_idx-piecepos);
    seg_idx = seg_idx>=0;
    seg_idx = sum(double(seg_idx));
    ax_coef = coef_struc.ax{seg_idx,1};
    ay_coef = coef_struc.ay{seg_idx,1};
    bx_coef = coef_struc.bx{seg_idx,1};
    by_coef = coef_struc.by{seg_idx,1};
    diff_z_mean = mean(coef_struc.diff_z(i,:));
    raw_data.z(idx) = raw_data.z(idx)-diff_z_mean; % z_coordinates are registered to m0 sequence
    % correct current z to the origin defined by the head of hyperstack segment, to which the current stack belongs 
    supp_z_ori=drift_z_m0(1+(piecepos(seg_idx)-1)*vol_per_hyper,:);
        
    drift_z_cur=drift_z_m0(i,:)'-supp_z_ori';
    
    [supp_z_ori,rank_order]=sort(supp_z_ori);
    
    drift_z_cur = drift_z_cur(rank_order);
    
    coef_z = polyfit(supp_z_ori,drift_z_cur,2);
    a = coef_z(1);
    b = coef_z(2);
    c = coef_z(3);
    
    z_obs = raw_data.z(idx);% z_obs-z_drift=z_cor    z_drift=a*z_cor^2+b*z_cor+c  z_cor is z support_ori
    z_0 = (-(b+1)+sqrt((b+1)^2-4*a*(c-z_obs)))/2/a;
    z_1 = (-(b+1)-sqrt((b+1)^2-4*a*(c-z_obs)))/2/a;
    dist_0 = abs(z_0-z_obs);
    dist_1 = abs(z_1-z_obs);
    if a~=0
        z_cor = zeros(size(z_obs));
        z_cor(dist_0<dist_1) = z_0(dist_0<dist_1);
        z_cor(dist_0>dist_1) = z_1(dist_0>dist_1);
    else
        z_cor = z_obs;
    end
    
    % z_cor is the corrected single molecule axial position in specified
    % hyperstack coordinate
    a_x = ax_coef(1)*z_cor + ax_coef(2);
    a_y = ay_coef(1)*z_cor + ay_coef(2);
    b_x = bx_coef(1)*z_cor.^2+bx_coef(2)*z_cor+bx_coef(3);
    b_y = by_coef(1)*z_cor.^2+by_coef(2)*z_cor+by_coef(3);
    
    t_cur = i-(piecepos(seg_idx)-1)*vol_per_hyper;
    delta_x = a_x*t_cur + b_x;
    delta_y = a_y*t_cur + b_y;
    raw_data.x(idx) = raw_data.x(idx)-delta_x;
    raw_data.y(idx) = raw_data.y(idx)-delta_y;
end
registered_data = raw_data;
end