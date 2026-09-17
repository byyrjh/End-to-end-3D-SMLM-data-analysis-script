function out = ls_mech_drift_FM(in)
vol_num = size(in.LS_os,1);
mat_coeff = zeros(vol_num,3);
% model ax+by+c = z
for i=1:vol_num
    LS_plane_data = in.LS_os(i,:);
    x_bin_pos = in.x(i,:)/100/10; % 100 is pixel size 10 is bin size
    y_bin_pos = in.y(i,:)/100/10; % 100 is pixel size 10 is bin size
    A = [x_bin_pos-1; y_bin_pos-1; ones(size(y_bin_pos))]';
    temp_coeff = inv(A'*A)*(A')*(LS_plane_data');
    mat_coeff(i,:) = temp_coeff';
end
for i = 1:3
    mat_coeff(:,i) = smooth(mat_coeff(:,i),50);
end
trace_os_m0_fitted = [];
for i=1:size(in.LS_os,2)
    x_bin_pos = in.x(:,i)/100/10-1; % 100 is pixel size 10 is bin size
    y_bin_pos = in.y(:,i)/100/10-1; % 100 is pixel size 10 is bin size
    z_fitted = x_bin_pos.*mat_coeff(:,1) + y_bin_pos.*mat_coeff(:,2) + mat_coeff(:,3);
    trace_os_m0_fitted = [trace_os_m0_fitted z_fitted-z_fitted(1)];
end
out = trace_os_m0_fitted;
end 