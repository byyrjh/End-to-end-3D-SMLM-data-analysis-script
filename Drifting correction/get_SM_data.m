function loc_SM = get_SM_data(dir, scan, setting_para,coef)
pixel_size_z = setting_para.pixel_size_z;  % nm
pixel_size_xy = setting_para.pixel_size_xy;  % nm
step_size = setting_para.step_size;  % nm
num_photon_factor = setting_para.num_photon_factor;
xybinsize_SM = setting_para.xybinsize; % unit nm
load(strcat(dir,'fitting_result_SM_',scan,'_round2.mat'));
load(strcat(dir,'map_ptr_t_SM_',scan,'.mat'));
load(strcat(dir,'map_ptr_x_SM_',scan,'.mat'));
load(strcat(dir,'map_ptr_y_SM_',scan,'.mat'));
load(strcat(dir,'map_ptr_z_SM_',scan,'.mat'));
% load(strcat(dir,'SM_os_aver_',scan,'.mat'))  % old version time axis smoothed light sheet offset map
% load(strcat(dir,'SM_os_AO_map_',scan,'.mat')) % old version Gaussian filter smoothed light sheet offset map following time trace smoothing
x_outlier_idx = (fitting_results(1,:)>4)';% single molecule is found at the edge of segment
y_outlier_idx = (fitting_results(2,:)>4)';% single molecule is found at the edge of segment
z_outlier_idx = (fitting_results(3,:)>40)';% single molecule is found far away from center slice of segment
os_outlier_idx = (abs(fitting_results(6,:))>60)';% light sheet offset is too large state transit possibly occurred in the segment
x_pos_SM=double((map_ptr_x_SM+fitting_results(1,:)')*pixel_size_xy);
y_pos_SM=double((map_ptr_y_SM+fitting_results(2,:)')*pixel_size_xy);
% assumption 1. Focal plane change resulting from RI fluctuation is
% negligible for dark red channnel. Therefore one time 5 parameters fitting
% is sufficient.
% assumption 2. Light sheet was stable during measurement that misalignment
% did not occur 
% assumption 3. Dynamic light sheet misalignment was considered in the way
% that two constant light sheet offsets were used for one time m0 and m1
% fitting with os_m1-os_m0 = dynamic misalignment
% assumption 4. Since static light sheet misalignment was unknown, |os_m1|
% == |os_m0| was assumed thus there was a undetermined yet globally constant
% localization bias for two channels
z_pos_SM=double(map_ptr_z_SM*step_size+(fitting_results(3,:)')*pixel_size_z); 
ls_os = double(fitting_results(6,:)')*pixel_size_z;
inten_SM=double(fitting_results(4,:)'*num_photon_factor);
t_SM=double(map_ptr_t_SM);
nan_x=isnan(x_pos_SM);
nan_y=isnan(y_pos_SM);
nan_z=isnan(z_pos_SM);
nan_os = isnan(ls_os);
kick_out_idx=nan_x|nan_y|nan_z|x_outlier_idx|y_outlier_idx|z_outlier_idx|os_outlier_idx|nan_os;
x_pos_SM(kick_out_idx)=[];
y_pos_SM(kick_out_idx)=[];
z_pos_SM(kick_out_idx)=[];
inten_SM(kick_out_idx)=[];
ls_os(kick_out_idx)=[];
t_SM(kick_out_idx)=[];

%%%%%%%%%%%%%  following linear model should be rechecked after FM channel
%%%%%%%%%%%%%  being corrected in the same way
%{
% calculate RI induced aberration
[xy t_range] = size(LS_os_AO_map);
xy = sqrt(xy);
LS_os_map = zeros(xy,xy,t_range);
LS_correction_map = zeros(xy,xy,t_range);
for i = 1:t_range
    LS_os_map(:,:,i) = reshape(LS_os_AO_map(:,i),[xy xy]);
end
mat_coeff = zeros(t_range,3);
aver_fit_err = zeros(t_range,1);
length_fit_data = zeros(t_range,1);
fit_data = cell(t_range,1);
% model ax+by+c = z
for i=1:t_range
    LS_plane_data = LS_os_map(:,:,i);
    x = [];
    y = [];
    z = [];
    for j=1:xy
        for k=1:xy
            if ~isnan(LS_plane_data(k,j))
                x = [x k];
                y = [y j];
                z = [z LS_plane_data(k,j)];
            end
        end
    end
    fit_data{i,1}.x = x;
    fit_data{i,1}.y = y;
    fit_data{i,1}.z = z;
    A = [x; y; ones(size(x))]';
    temp_coeff = inv(A'*A)*(A')*(z');
    mat_coeff(i,:) = temp_coeff';
    z_fit = (x)*temp_coeff(1)+(y)*temp_coeff(2)+temp_coeff(3);
    aver_fit_err(i) = sqrt(sum((z_fit-z).^2)/length(x));
    length_fit_data(i) = length(x);
end

for y=1:xy
    for x=1:xy
        data_raw_temp = squeeze(LS_os_map(x,y,:));
        data_fit = (x)*mat_coeff(:,1)+(y)*mat_coeff(:,2)+mat_coeff(:,3);
        os_diff = data_raw_temp-data_fit;
        LS_correction_map(x,y,:) = reshape(os_diff,[1 1 t_range]);
    end
end
%}

if strcmp(scan,'m1') 
    % z_m1_SM registration -> z_m0_reg_SM(t) = z_m1_apparent_SM(t) + ls_os_m1_apparent_SM(t) - z_diff(t) 
    z_m0_reg_SM = zeros(size(z_pos_SM));
    x_m0_reg_SM = zeros(size(x_pos_SM));
    y_m0_reg_SM = zeros(size(y_pos_SM));
    for i = 1:length(coef.diff_z)
        t_idx = (t_SM == i);
        % register and pre-correct z 
        z_m0_reg_SM(t_idx) = z_pos_SM(t_idx) - coef.diff_z(i) + ls_os(t_idx);
        % register x and y
        diff_x = z_m0_reg_SM(t_idx).^2*coef.x(i,1) + z_m0_reg_SM(t_idx)*coef.x(i,2) + coef.x(i,3); 
        diff_y = z_m0_reg_SM(t_idx).^2*coef.y(i,1) + z_m0_reg_SM(t_idx)*coef.y(i,2) + coef.y(i,3); 
        x_m0_reg_SM(t_idx) = x_pos_SM(t_idx) - diff_x;
        y_m0_reg_SM(t_idx) = y_pos_SM(t_idx) - diff_y;
    end 
    loc_SM.x=x_m0_reg_SM;
    loc_SM.y=y_m0_reg_SM;
    loc_SM.z=z_m0_reg_SM;
    loc_SM.inten=inten_SM;
    loc_SM.t=t_SM;

    %{
    load(strcat(dir,'SM_os_AO_4D_map_m0.mat'));
    SlicesPerVol = setting_para.SlicesPerVol;
    [xyz t] = size(LS_os_AO_4D_map);
    xy = sqrt(xyz/SlicesPerVol); 
    for i = 1:length(z_pos_SM)
        cur_x_bin = round(x_pos_SM(i) / pixel_size_xy) / (xybinsize_SM / pixel_size_xy);
        cur_x_bin = ceil(cur_x_bin);
        cur_y_bin = round(y_pos_SM(i) / pixel_size_xy) / (xybinsize_SM / pixel_size_xy);
        cur_y_bin = ceil(cur_y_bin);
        cur_z = ceil(z_pos_SM(i)/step_size);
        cur_t = t_SM(i);
        idx_sp = (cur_z-1)*xy*xy + (cur_y_bin-1)*xy + cur_x_bin;
        temp_ls_os_m0 = LS_os_AO_4D_map(idx_sp,cur_t);
        if (~isnan(temp_ls_os_m0))
% the logic is similar to diff_z 
% diff_z = drift_z_m1 - drift_z_m0 + (FM_trace_m1.LS_os-FM_trace_m0.LS_os);
% drift_z_m0 = drift_z_m1 + (FM_trace_m1.LS_os-FM_trace_m0.LS_os) - diff_z
% the idea here is to convert apparent loc_z_m1 to apparent loc_z_m0
% followed by pre_correct_z function that converts apparent loc_z_m0 to
% physical loc_z
% in reality -temp_ls_os_m0 will be cancelled off in pre_correct_z function
            z_pos_SM(i) = z_pos_SM(i) + (ls_os(i)-temp_ls_os_m0)*setting_para.pixel_size_z;
        end
    end
    %}
else 
    % pre-correct z
    % z_m0_SM = z_m0_apparent_SM + ls_os_m0_apparent_SM
    z_m0_pre_cor = z_pos_SM + ls_os;
    loc_SM.x=x_pos_SM;
    loc_SM.y=y_pos_SM;
    loc_SM.z=z_m0_pre_cor;
    loc_SM.inten=inten_SM;
    loc_SM.t=t_SM;
end

end