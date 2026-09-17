function [para_quadratic_fit FM_trace_m0]= get_correction_coef(dir,setting_para)
smooth_f = setting_para.smooth_f;% localization precision is independent of smooth binning size
exposuretime = setting_para.exposuretime;  %ms
slice_per_stack = setting_para.slice_per_stack;
coef_vis = setting_para.coef_vis;
diff_vis = setting_para.diff_vis;
scan_m0 = 'm0';
scan_m1 = 'm1';
FM_trace_m0 = get_FM_trace(dir, scan_m0, setting_para);
FM_trace_m1 = get_FM_trace(dir, scan_m1, setting_para);

%% debug
figure
hold on
set(gca,'FontSize',20)
scatter3(FM_trace_m0.x(1,:)/1000,FM_trace_m0.y(1,:)/1000,FM_trace_m0.z(1,:)/1000,'LineWidth',4);%,'MarkerEdgeColor',color_map_list(i,:));
scatter3(FM_trace_m1.x(1,:)/1000,FM_trace_m1.y(1,:)/1000,FM_trace_m1.z(1,:)/1000,'LineWidth',4);%,'MarkerEdgeColor',color_map_list(i,:));
axis equal
xlabel('scanning / m')
ylabel('propagation / m')
zlabel('z / m')
box on
grid on
title('Fiducial marker distribution')
view(90,0);
%%
mkdir(strcat(dir,'register_vis'))
[FM_trace_m0 FM_trace_m1] = FM_pair(FM_trace_m0, FM_trace_m1, setting_para);
[t_length, num_FM_tr] = size(FM_trace_m0.x);
vol_per_hyper = setting_para.vol_per_hyper-setting_para.stationary_pos;
num_hyperstack = setting_para.num_hyper;

drift_x_m0 = smooth_gen(FM_trace_m0.x, vol_per_hyper, setting_para.piecepos, smooth_f);
drift_y_m0 = smooth_gen(FM_trace_m0.y, vol_per_hyper, setting_para.piecepos, smooth_f);
drift_z_m0 = smooth_gen(FM_trace_m0.z, vol_per_hyper, setting_para.piecepos, smooth_f);
drift_x_m1 = smooth_gen(FM_trace_m1.x, vol_per_hyper, setting_para.piecepos, smooth_f);
drift_y_m1 = smooth_gen(FM_trace_m1.y, vol_per_hyper, setting_para.piecepos, smooth_f);
drift_z_m1 = smooth_gen(FM_trace_m1.z, vol_per_hyper, setting_para.piecepos, smooth_f);

LS_os_rel_m0 = FM_trace_m0.LS_os-FM_trace_m0.LS_os(1,:);
LS_os_rel_m1 = FM_trace_m1.LS_os-FM_trace_m1.LS_os(1,:);
diff_x = drift_x_m1-drift_x_m0;
diff_y = drift_y_m1-drift_y_m0;
% diff_z = drift_z_m1-(drift_z_m0);
%%%%%%%%%%%%% the following correction method is only valid provided the
%%%%%%%%%%%%% lightsheet is free of mechanical misalignement, which is very
%%%%%%%%%%%%% rare given due to the poor stability of mirror mounts and SLM
%%%%%%%%%%%%% mount, that gives rise to the hight RMSE compared to the
%%%%%%%%%%%%% result without correct (see running results of FM_loc_prec.m)
% diff_z = drift_z_m1+LS_os_rel_m1-(drift_z_m0+LS_os_rel_m0); % mechanical diff is independent of wavelength 

%%%%%%%%%%%%%%%%%%%   11.1.2024   %%%%%%%%%%%%%%%%%%%%
diff_z = drift_z_m1 - drift_z_m0 + (FM_trace_m1.LS_os-FM_trace_m0.LS_os);
% diff_light sheet offset is considered to be constituted from
% aberration induced offset and dynamics misalignment (phase delay). Both
% components result in detective objective focal plane change, therefore
% needed to be subtracted from localization in z
% the break of covarying of diff_ls_os is the consequence of aberration,
% which distorts the light sheet plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%    5.3.2024    %%%%%%%%%%%%%%%%%%%
% light sheet misalignment ax + by + c = z 
% for m0 and m1 scan, a_m0 b_m0 should equal to a_m1 and b_m1
% c_m0 differs from c_m1 by constant dynamic misalignment
% model ax+by+c = z  a marks scanning direction and b marks propogation
%        if light sheet offset is globally changing
%        one cannot tell if this change is the global misalignment or
%        global RI change induced focal plane change. Whereas the RI induced 
%        light sheet offset change should be added to z localization, the
%        mechanical misalignment shouldn't contribute to z localization
%        errorr
%        RI change not only causes irregular fluctuation of ls os, but also
%        contributes to glolal light sheet plane tip-and-tilt and shifting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this Diff is deemed to be global constant consisted of hysteresis and phase delay due to high frequency scan 
% the reason why it appears to be a tilt line is the consequence of
% nonlinear hysteresis and z drifting
drift_z_m0 = drift_z_m0;% +LS_os_rel_m0; drift_z_m0 = drift_z_m0 +LS_os_rel_m0; this causes unpredicted time trace
% drift_z_m0 = drift_z_m0;
% diff_z = drift_z_m1+FM_trace_m1.LS_os-(drift_z_m0+FM_trace_m0.LS_os);
% Owing to dynamic misalignment the apparent z_m0 localization differs from z_m1 by distance equaling to dynamic misalignment, 
% which also applies to dark red channel provided correct light sheet
% offsets are assigned to one round 5 parameter fitting in dark red channel
diff_os = FM_trace_m1.LS_os-FM_trace_m0.LS_os;
% this constant Diff independent of t speaks of a static dynamic
% misalignment (constant phase delay)
% Because of the existance of hysteresis in z direction diff_os varies with
% depth where FMs are found despite the constant lag shown on oscilloscope,
% meaning the position sensor on piezo stage is not sufficiently sensitive


% LS_os_m0 = single(mean(FM_trace_m0.LS_os,2)/10);
% LS_os_m1 = single(mean(FM_trace_m1.LS_os,2)/10);
% save(strcat(dir,'LS_os_m0.mat'),'LS_os_m0');
% save(strcat(dir,'LS_os_m1.mat'),'LS_os_m1');

if diff_vis
    color_map_list = cm_gen(num_FM_tr);
    [~, z_rank_plot] = sort(drift_z_m0(1,:)/1000);
    vol_time = exposuretime*slice_per_stack*2;
    tot_vol = vol_per_hyper*num_hyperstack;
    t_axis = linspace(0,tot_vol*vol_time-vol_time,tot_vol)/1000/60;
    hyper_pos = linspace(0,tot_vol*vol_time-vol_time*vol_per_hyper,num_hyperstack)/1000/60;
    % light sheet offset
    figure
    hold on
    grid on
    % title('Light sheet offset')
    ylabel('Distance / nm')
    xlabel('t / min')
    LS_os_vis = FM_trace_m0.LS_os;
    LS_os_vis = LS_os_vis(:,z_rank_plot);
    for i=1:num_FM_tr
        plot(t_axis,LS_os_vis(:,i),'Color',color_map_list(i,:));
%         plot(t_axis,FM_trace_m1.LS_os(:,i))
    end
    set(gca,'FontSize',18)
    saveas(gcf,strcat(dir,'register_vis\','Light sheet offset.bmp'))
    % difference x
    diff_x_vis = diff_x(:,z_rank_plot);
    figure
    hold on
    for i=1:num_FM_tr
        plot(t_axis,(diff_x_vis(:,i)),'-','Color',color_map_list(i,:))
    end
    for i=1:num_hyperstack
        if mod(i,10) == 0
            xline(hyper_pos(i),'LineWidth',1)
        else
            xline(hyper_pos(i),'--','LineWidth',0.5)
        end
    end
    % title('Difference in x')
    set(gca,'FontSize',20)
    grid off
    box on
    xlabel('t / min')
    ylabel('Distance / nm')
    saveas(gcf,strcat(dir,'register_vis\','Difference in x.bmp'))

    % difference y
    diff_y_vis = diff_y(:,z_rank_plot);
    figure
    hold on
    for i=1:num_FM_tr
        plot(t_axis,(diff_y_vis(:,i)),'-','Color',color_map_list(i,:))
    end
    for i=1:num_hyperstack
        if mod(i,10) == 0
            xline(hyper_pos(i),'LineWidth',1)
        else
            xline(hyper_pos(i),'--','LineWidth',0.5)
        end
    end
    % title('Difference in propagation direction')
    set(gca,'FontSize',20)
    grid off
    box on
    xlabel('t / min')
    ylabel('Distance / nm')
    saveas(gcf,strcat(dir,'register_vis\','Difference in propagation direction.bmp'))
%%                   2025_5_22      implemented 
%%%%%%%%%%%%%%%   Andrei workflow   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = (drift_x_m1(1,:) + drift_x_m0(1,:))/2;
y0 = (drift_y_m1(1,:) + drift_y_m0(1,:))/2;
for i = 1:8
drift_x_m1_orig0(:,i) = drift_x_m1(:,i) - x0(i);
drift_x_m0_orig0(:,i) = drift_x_m0(:,i) - x0(i);
drift_y_m1_orig0(:,i) = drift_y_m1(:,i) - y0(i);
drift_y_m0_orig0(:,i) = drift_y_m0(:,i) - y0(i);
end
x_mech = [drift_x_m0_orig0 drift_x_m1_orig0];
x_mech = mean(x_mech,2);

figure
hold on 
for i = 1:8
    plot(drift_x_m0_orig0(:,i),'Color',color_map_list(i,:));
    plot(drift_x_m1_orig0(:,i),'Color',color_map_list(i,:));
end
plot(x_mech,'-k',LineWidth=2)

Res_x_m0 = drift_x_m0_orig0-x_mech;
Res_x_m1 = drift_x_m1_orig0-x_mech;
Res_x_m0_mean = mean(Res_x_m0,2);
Res_x_m1_mean = mean(Res_x_m1,2);
figure
hold on 
for i = 1:8
    plot(Res_x_m0(:,i),'Color',color_map_list(i,:));
    plot(Res_x_m1(:,i),'Color',color_map_list(i,:));
end
plot(Res_x_m0_mean,'-k',LineWidth=2)
plot(Res_x_m1_mean,'-k',LineWidth=2)


z0 = (drift_z_m1(1,:) + drift_z_m0(1,:))/2;
for i = 1:8
drift_z_m1_orig0(:,i) = drift_z_m1(:,i) - z0(i);
drift_z_m0_orig0(:,i) = drift_z_m0(:,i) - z0(i);
end
z_mech = [drift_z_m0_orig0 drift_z_m1_orig0];
z_mech = mean(z_mech,2);
figure
hold on 
for i = 1:8
    plot(drift_z_m0_orig0(:,i),'Color',color_map_list(i,:));
    plot(drift_z_m1_orig0(:,i),'Color',color_map_list(i,:));
end
plot(z_mech,'-k',LineWidth=2)


%%%%%%%%%%%%%%%   Andrei workflow   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    % difference z
    diff_z_vis = diff_z(:,z_rank_plot);
    figure
    hold on
    for i=1:num_FM_tr
        plot(t_axis,smooth(diff_z_vis(:,i),30),'-','Color',color_map_list(i,:))
        diff_z_vis_save(:,i) = smooth(diff_z_vis(:,i),30);
    end
    plot(t_axis,mean(diff_z,2),'LineWidth',2,'Color','k');
    for i=1:num_hyperstack
        if mod(i,10) == 0
            xline(hyper_pos(i),'LineWidth',1)
        else
            xline(hyper_pos(i),'--','LineWidth',0.5)
        end
    end
    % title('Difference in z')
    set(gca,'FontSize',20)
    grid off
    box on
    xlabel('t / min')
    ylabel('$\textrm{z}_\textrm{mb}+\bar{\textrm{z}}_\textrm{mb}-(\textrm{z}_\textrm{mf}+\bar{\textrm{z}}_\textrm{mf})$ / nm','Interpreter','latex')
%     ylim([-400, -200])
    saveas(gcf,strcat(dir,'register_vis\','Difference in z.bmp'))

    % for comp
    figure
    diff_z_wo_ls_correction = drift_z_m1 - drift_z_m0;
    diff_z_wo_ls_correction_vis = diff_z_wo_ls_correction(:,z_rank_plot);
    hold on
    for i=1:num_FM_tr
        plot(t_axis,smooth(diff_z_wo_ls_correction_vis(:,i),30),'-','Color',color_map_list(i,:))
        diff_z_wo_ls_correction_vis_save(:,i) = smooth(diff_z_wo_ls_correction_vis(:,i),30);
    end
    plot(t_axis,mean(diff_z_wo_ls_correction,2),'LineWidth',2,'Color','k');
    for i=1:num_hyperstack
        if mod(i,10) == 0
            xline(hyper_pos(i),'LineWidth',1)
        else
            xline(hyper_pos(i),'--','LineWidth',0.5)
        end
    end
    % title('Difference in z without ls correction')
    set(gca,'FontSize',20)
    grid off
    box on
    xlabel('t / min')
    ylabel('$\textrm{z}_\textrm{mb}-\textrm{z}_\textrm{mf}$ / nm','Interpreter','latex')
    saveas(gcf,strcat(dir,'register_vis\','Difference in z without ls correction.bmp'))
    
    % absolute light sheet offset m0
    figure
    hold on
    for i=1:num_FM_tr
        plot(t_axis,(LS_os_vis(:,i)),'-','Color',color_map_list(i,:))
    end
    plot(t_axis,mean(FM_trace_m0.LS_os,2),'LineWidth',2,'Color','k')
    for i=1:num_hyperstack
        if mod(i,10) == 0
            xline(hyper_pos(i),'LineWidth',1)
        else
            xline(hyper_pos(i),'--','LineWidth',0.5)
        end
    end
    % title('Absolute LS offset m0')
    set(gca,'FontSize',20)
    grid off
    box on
    xlabel('t / min')
    ylabel('$\bar{\textrm{z}}_\textrm{mf}$ / nm','Interpreter','latex')

    % difference light sheet offset
    diff_os_vis = diff_os(:,z_rank_plot);
    figure
    hold on
    for i=1:num_FM_tr
        plot(t_axis,(diff_os_vis(:,i)),'-','Color',color_map_list(i,:))
    end
    plot(t_axis,mean(diff_os,2),'LineWidth',2,'Color','k')
    for i=1:num_hyperstack
        if mod(i,10) == 0
            xline(hyper_pos(i),'LineWidth',1)
        else
            xline(hyper_pos(i),'--','LineWidth',0.5)
        end
    end
    % title('Difference in LS offset')
    set(gca,'FontSize',20)
    grid off
    box on
    xlabel('t / min')
    ylabel('$\bar{\textrm{z}}_\textrm{mb}-\bar{\textrm{z}}_\textrm{mf}$ / nm','Interpreter','latex')
    saveas(gcf,strcat(dir,'register_vis\','Difference in LS_os.bmp'))
end

drift_x_m0_vis = drift_x_m0(1,z_rank_plot)/1000;
drift_y_m0_vis = drift_y_m0(1,z_rank_plot)/1000;
drift_z_m0_vis = drift_z_m0(1,z_rank_plot)/1000;

figure
hold on
set(gca,'FontSize',20)
for i=1:size(drift_x_m0,2)
scatter3(drift_x_m0_vis(i),drift_y_m0_vis(i),drift_z_m0_vis(i),'LineWidth',4,'MarkerEdgeColor',color_map_list(i,:));
end
axis equal
xlabel('scanning / m')
ylabel('propagation / m')
zlabel('z / m')
box on
grid on
title('Fiducial marker distribution')
view(90,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     diff_x and diff_y are piecewise funtions      %%%%%%%%%%%%%
%%%     diff_x_1 = a*t + b
%%%     diff_x_1 = (coef_ax_1(1)*z + coef_ax_1(2))*t+(coef_bx_1(1)*z^2+coef_bx_1(2)*z+coef_bx_1(3))
%%%     z is the initial value at the head of each FM_m0 piece
%%%%%%%%%%%%%%%%%%%%%%%%%%%  6.3.2024  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     z_m0_physical_FM = z_m0_apparent_FM + ls_os_m0_apparent_FM 
%%%     z_m1_physical_FM = z_m1_apparent_FM + ls_os_m1_apparent_FM 
%%%     z_diff = z_m1_physical_FM - z_m0_physical_FM
%%%     z_m1_SM registration -> z_m0_reg_SM = z_m1_apparent_SM + ls_os_m1_apparent_SM - z_diff
%%%     z_m0_SM = z_m0_apparent_SM + ls_os_m0_apparent_SM
%%%     z_m0_reg_SM and z_m0_SM are treated equally in drifting correction for z
%%%     drifting correction: z_m0_SM_real(t) = z_m0_SM(t) - (z_m0_physical_FM(t) - z_m0_physical_FM(0))
%%
para_quadratic_fit.x = zeros(size(diff_x,1),3);
para_quadratic_fit.y = zeros(size(diff_y,1),3);
for i = 1:size(diff_x,1)
    % fitting support coordinate is z_m0_physical_FM(t)
    support_coordinate_z = drift_z_m0(i,:) + FM_trace_m0.LS_os(i,:);
    para_quadratic_fit.x(i,:) = polyfit(support_coordinate_z,diff_x(i,:),2);
    para_quadratic_fit.y(i,:) = polyfit(support_coordinate_z,diff_y(i,:),2);
end
para_quadratic_fit.diff_z = mean(diff_z,2);
%% manual piecewise fitting
% see old version
end