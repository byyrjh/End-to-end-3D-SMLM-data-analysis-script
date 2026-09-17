function [plane_fit_para plane_fit_map] = ls_plan_fit(in)
os_unit = 10; % nm
later_unit = 1; % um
smooth_f_z = 3;
smooth_f_t = 100;
[x y z t] = size(in);
mat_coeff = zeros(3,z,t);
xx = -(x-1)/2:later_unit:(x-1)/2;
[xx yy]= meshgrid(xx,xx);
in = in * os_unit;
parfor k = 1:t
    xx_temp = xx;
    yy_temp = yy;
    cur_vol_data = in(:,:,:,k);
    for j = 1:z
        cur_plane_data = cur_vol_data(:,:,j);
        idx = ~isnan(cur_plane_data);
        cur_plane_data = cur_plane_data(idx);
        temp_x = xx_temp(idx);
        temp_y = yy_temp(idx);
        A = [temp_x temp_y ones(size(temp_y))];
        try
        temp_coeff = (A'*A)\(A')*(cur_plane_data);
        mat_coeff(:,j,k) = temp_coeff;
        catch
        end
    end
end
for i = 1:z
    for j = 1:3
        temp = squeeze(mat_coeff(j,i,:));
        temp = smooth(temp,smooth_f_t);
        mat_coeff(j,i,:) = temp;
    end
end
for i = 1:t
    for j = 1:3
        temp = squeeze(mat_coeff(j,:,i));
        temp = smooth(temp,smooth_f_z);
        mat_coeff(j,:,i) = temp;
    end
end
plane_fit_para = mat_coeff;
plane_fit_map = zeros(size(in));
for i = 1:t
    for j = 1:z
        temp_a = mat_coeff(1,j,i);
        temp_b = mat_coeff(2,j,i);
        temp_c = mat_coeff(3,j,i);
        plane_fit_map(:,:,j,i) = xx*temp_a+yy*temp_b+temp_c;
    end
end
plane_fit_map = plane_fit_map-plane_fit_map(:,:,:,1);
end