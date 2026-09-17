function [FM_trace_m0_new FM_trace_m1_new] = FM_pair(FM_trace_m0, FM_trace_m1, setting_para)
%%
% FM trajectories are paired and sorted in ascending order
origin_drifting = setting_para.origin_drifting;
ori_m0 = [FM_trace_m0.x(origin_drifting,:)' FM_trace_m0.y(origin_drifting,:)' FM_trace_m0.z(origin_drifting,:)'];
ori_m1 = [FM_trace_m1.x(origin_drifting,:)' FM_trace_m1.y(origin_drifting,:)' FM_trace_m1.z(origin_drifting,:)'];
[num_FM_m0,~] = size(ori_m0);
[num_FM_m1,~] = size(ori_m1);
reoccur_m0 = [];
reoccur_m1 = [];
err_flag = false;
for i = 1:num_FM_m0
    for j = 1:num_FM_m1
        dist = sqrt(sum((ori_m0(i,:)-ori_m1(j,:)).^2));
        if (dist<setting_para.backlash_thres)
            err_flag = err_flag | sum(reoccur_m0==i)>0;
            err_flag = err_flag | sum(reoccur_m1==j)>0;
            reoccur_m0 = [reoccur_m0 i];
            reoccur_m1 = [reoccur_m1 j];
        end
    end
end
if err_flag
    error('Please set smaller backlash distance')
end
z = FM_trace_m0.z(origin_drifting,reoccur_m0);
[z_rank idx] = sort(z);
reoccur_m0 = reoccur_m0(idx);
reoccur_m1 = reoccur_m1(idx);
FM_trace_m0_new.x = FM_trace_m0.x(:,reoccur_m0);
FM_trace_m0_new.y = FM_trace_m0.y(:,reoccur_m0);
FM_trace_m0_new.z = FM_trace_m0.z(:,reoccur_m0);
FM_trace_m0_new.inten = FM_trace_m0.inten(:,reoccur_m0);
FM_trace_m0_new.LS_os = FM_trace_m0.LS_os(:,reoccur_m0);

FM_trace_m1_new.x = FM_trace_m1.x(:,reoccur_m1);
FM_trace_m1_new.y = FM_trace_m1.y(:,reoccur_m1);
FM_trace_m1_new.z = FM_trace_m1.z(:,reoccur_m1);
FM_trace_m1_new.inten = FM_trace_m1.inten(:,reoccur_m1);
FM_trace_m1_new.LS_os = FM_trace_m1.LS_os(:,reoccur_m1);
end