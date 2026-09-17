function corrected_data = correct_m0(raw_data,raw_drifting_tr,setting_para)
axis_z_ori=raw_drifting_tr.z(1,:);
[axis_z_ori,rank_order]=sort(axis_z_ori);
t_SM = raw_data.t;
num_vol = max(raw_data.t);
raw_drifting_tr_x = raw_drifting_tr.x;
raw_drifting_tr_y = raw_drifting_tr.y;
raw_drifting_tr_z = raw_drifting_tr.z;
LS_os_rel = raw_drifting_tr.LS_os_rel;
raw_data = pre_correct_z(raw_data,setting_para);
for i=1:num_vol
        idx=(t_SM==i);
        drift_x_cur=raw_drifting_tr_x(i,:)'-raw_drifting_tr_x(1,:)';
        drift_y_cur=raw_drifting_tr_y(i,:)'-raw_drifting_tr_y(1,:)';
        drift_z_cur=raw_drifting_tr_z(i,:)'-raw_drifting_tr_z(1,:)';
        sup_z_cur = raw_drifting_tr_z(1,rank_order)';
        LS_os_cur = LS_os_rel(i,:)';
        drift_x_cur=drift_x_cur(rank_order);
        drift_y_cur=drift_y_cur(rank_order);
        drift_z_cur=drift_z_cur(rank_order)+LS_os_cur(rank_order);
        coef_x = polyfit(axis_z_ori,drift_x_cur,2);
        coef_y = polyfit(axis_z_ori,drift_y_cur,2);
        coef_z = polyfit(sup_z_cur,drift_z_cur,2);
        a = coef_z(1);
        b = coef_z(2);
        c = coef_z(3);
        % two method (1 initial positions as support set, solve quadratic equation 2 observed position as support set, direct substitute)
        % correct for light sheet offset induced bias(originated from focal position change)
        z_obs = raw_data.z(idx);% z_obs-z_drift=z_cor    z_drift=a*z_cor^2+b*z_cor+c
        z_0 = (-(b+1)+sqrt((b+1)^2-4*a*(c-z_obs)))/2/a;
        z_1 = (-(b+1)-sqrt((b+1)^2-4*a*(c-z_obs)))/2/a;
        dist_0 = abs(z_0-z_obs);
        dist_1 = abs(z_1-z_obs);
%         raw_data.z(idx) = cur_z_cor-(coef_z(1)*cur_z_cor.^2+coef_z(2)*cur_z_cor+coef_z(3));
        if i~=1
            z_cor = zeros(size(z_obs));
            z_cor(dist_0<dist_1) = z_0(dist_0<dist_1);
            z_cor(dist_0>dist_1) = z_1(dist_0>dist_1);
%             raw_data.z(idx) = z_cor;
        end
        raw_data.z(idx) = raw_data.z(idx) - mean(drift_z_cur);
        raw_data.x(idx)=raw_data.x(idx)-(coef_x(1)*raw_data.z(idx).^2+coef_x(2)*raw_data.z(idx)+coef_x(3));
        raw_data.y(idx)=raw_data.y(idx)-(coef_y(1)*raw_data.z(idx).^2+coef_y(2)*raw_data.z(idx)+coef_y(3));
end
corrected_data = raw_data;
end