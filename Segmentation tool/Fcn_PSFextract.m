function Fcn_PSFextract(hObj,event,obj)
% ADU_thres_diff=15;
% FM_flag=false;
%%
% parameter input volume half width, volume half height determine segment size
% parameter input half height in search determines half of search box depth
% paramter input prominence determines hierachy search position is prominence times of curent global peak 
% parameter input ADU_thres determines the ADU enforcing termination of searching
% parameter input std_filter determine stadard dceviation of Gaussian filter
% hardcode parameter Hbox_size determines half of search box width
% hardcode parameter ideal PSF lateral FWHM 4 pixels
% hardcode parameter ideal PSF axial FWHM is twice of lateral FWHM
% hardcode parameters 0.6 and 1 damp ADU down to 1/curpeak times at 0.6 pixel 
% radius in later and 1 pixel in axial area around localization spot
%%
% Laplacian filter reduces to 2D for low SNR cell data
% hard coded lateral search range 9*9 or 21*21
%%
Hbox_size_search_diff=0;
pararead=findobj('Tag','PSFpara');
width_in=findobj('Tag','width');
height_in=findobj('Tag','height');
height_search_in=findobj('Tag','height_search');
width_search_in=findobj('Tag','width_search');
ADU_thres_in=findobj('Tag','ADU_thres');
FM_label = findobj('Tag','FM_button_activation');
ana_all = findobj('Tag','cb_ana_all');
if FM_label.Value
    FM_flag=true;
    drift_ampli_in=findobj('Tag','drift_ampli_in');
    ADU_thres_diff = str2double(drift_ampli_in.String); 
    obj.guihandles.PSFdata.SegPara.datatype = 'FM';
    obj.guihandles.PSFdata.SegPara.ADU_diff = ADU_thres_diff;
else 
    obj.guihandles.PSFdata.SegPara.datatype = 'SM';
    FM_flag=false;
    ADU_thres_diff=0;
    obj.guihandles.PSFdata.SegPara.ADU_diff = ADU_thres_diff;
end

if ana_all.Value
    path_sample = obj.guihandles.PSFdata.path{1,1};
    slash_pos = strfind(path_sample,'\');
    path = path_sample(1:slash_pos(end));
    det = strfind(path_sample,'det');
    det = path_sample(det:det+3);
    scan = strfind(path_sample,'_m=');
    scan = path_sample(scan:scan+3);
    dir_info = dir(strcat(path,'*',scan,'*',det,'.tif'));
    nfile = size(dir_info,1);
    path_ary = cell(nfile,1);
    for i = 1 : nfile
        path_ary{i,1} = strcat(dir_info(i).folder,'\',dir_info(i).name);
    end
    obj.guihandles.PSFdata.path = path_ary;
else
    path_ary=obj.guihandles.PSFdata.path;
end


% Std_filter_in=findobj('Tag','Std_filter');

% parameters setting
width_seg=str2double(width_in.String);  % usually 30 with lattice light sheet microscopy configuration
height=str2double(height_search_in.String);  % unit in step
height_seg=str2double(height_in.String);  % unit in step
stepsize=obj.guihandles.PSFdata.stepsize;  % unit um
Hbox_size_xy=str2double(width_search_in.String); % search area is hard coded as 21*21 pixels
ADU_thres=str2double(ADU_thres_in.String); 
std = 1;%str2double(Std_filter_in.String);  % 3D Gaussian filter
pixelsize=(obj.guihandles.PSFdata.pixelsize)/(obj.guihandles.PSFdata.magnification)*obj.guihandles.PSFdata.binsize;  % unit um
ROI_sele=obj.guihandles.menu_analyze_ROIselect_flag.old;
obj.guihandles.PSFdata.slices=ROI_sele.zmax-ROI_sele.zmin+1;
binsize=obj.guihandles.PSFdata.binsize;
obj.guihandles.PSFdata.SegPara.width_in=width_seg;
obj.guihandles.PSFdata.SegPara.height_in=height_seg;
obj.guihandles.PSFdata.SegPara.width_search_in=Hbox_size_xy;
obj.guihandles.PSFdata.SegPara.height_search_in=height;
obj.guihandles.PSFdata.SegPara.ADU_thres_in=ADU_thres;
obj.guihandles.PSFdata.SegPara.Std_filter_in=std;
close(pararead);
%%   new algorithm parameter settings
pixel_size=pixelsize*1000;  % unit nm
step_size=stepsize*1000; % unit nm
filter_size=[5 5 5];
prominence=0.5; 
para.width_seg=width_seg;
para.height_seg=height_seg;
para.Hbox_size_axi=height; % unit in step 
para.Hbox_size_search=Hbox_size_xy;
para.FM=FM_flag;
Hbox_size_axi=height;
% global img3D_proc;
filter_size_in=(filter_size-1)/2;
[xx,yy,zz] = meshgrid(-filter_size_in(2):filter_size_in(2),-filter_size_in(1):filter_size_in(1),-filter_size_in(3):filter_size_in(3));
aspect_ratio=step_size/pixel_size;
filter_Gau=exp(-(xx.^2+yy.^2+(aspect_ratio*zz).^2)/(2*std^2));
filter_Gau(filter_Gau<eps*max(filter_Gau(:))) = 0;
para.filter_Gau=filter_Gau/sum(filter_Gau(:));
sigma=sqrt(2/log(2))/binsize;  %FWHM=4 pixels such that exp(-2^2/2/sigma^2)=1/2 without binning
sim_x=linspace(-Hbox_size_xy,Hbox_size_xy,Hbox_size_xy*2+1); % both sim_x and sim_z in unit of lateral pixel
sim_z=linspace(-Hbox_size_axi*aspect_ratio,Hbox_size_axi*aspect_ratio,Hbox_size_axi*2+1); 
[sim_xx sim_yy sim_zz]=meshgrid(sim_x,sim_x,sim_z);
para.sim_dist_std=0;
para.sim_median_dist=0;  % (median is independent of Hbox_size_axi and step size as simulated PSF is isotropic)???
para.sim_mean_dist=0;
sta_temp=[];
Gx_oprt=[-1 0 1];
Gy_oprt=[-1;0;1];
Gz_oprt=[-1;0;1];
%%
sim_ADU_center=exp(-(sim_xx.^2+sim_yy.^2+((sim_zz).^2)/4)/2/sigma^2);
sim_ADU_median_center=zeros(Hbox_size_xy*2+1,Hbox_size_xy*2+1,Hbox_size_axi*2+1);
    for k=1:Hbox_size_axi*2+1
        sim_ADU_median_center(:,:,k)=medfilt2(sim_ADU_center(:,:,k));
    end
    sim_ADU_median_Gau_center=imfilter(sim_ADU_median_center,para.filter_Gau);
    sim_ADU_median_Gau_Lap_center=div_img(sim_ADU_median_Gau_center,Gx_oprt,Gy_oprt,Gz_oprt,step_size);
    sim_ADU_center=sim_ADU_median_Gau_Lap_center;
    normal_factor=max(sim_ADU_center,[],'all');
%%
for i=1:step_size/5
    z_offset=(-step_size/5/2+i)*5/pixel_size;
    sim_ADU=exp(-(sim_xx.^2+sim_yy.^2+((sim_zz+z_offset).^2)/4)/2/sigma^2); % assuming PSF aspect ratio of 2
    sim_ADU_median=zeros(size(sim_ADU));
    for k=1:Hbox_size_axi*2+1
        sim_ADU_median(:,:,k)=medfilt2(sim_ADU(:,:,k));
    end
    sim_ADU_median_Gau=imfilter(sim_ADU_median,para.filter_Gau);
    sim_ADU_median_Gau_Lap=div_img(sim_ADU_median_Gau,Gx_oprt,Gy_oprt,Gz_oprt,step_size)/normal_factor;
    sim_ADU_raw=sim_ADU;
    sim_ADU=sim_ADU_median_Gau_Lap;       
%     sim_ADU(sim_ADU==1)=0;
    sim_idx=find(sim_ADU>prominence);
    sim_rr=sqrt(sim_xx.^2+sim_yy.^2+sim_zz.^2);
    sim_dist_arr=sim_rr(sim_idx);
    sim_dist_std=sqrt(sum(sim_dist_arr.^2)/length(sim_dist_arr));
    sta_temp=[sta_temp,length(sim_dist_arr)];
end
para.num_px_thres=min(sta_temp);

% para.sim_dist_std=sim_dist_std;
% para.sim_median_dist=median(sim_dist_arr);  % (median is independent of Hbox_size_axi and step size as simulated PSF is isotropic)???
% para.sim_mean_dist=mean(sim_dist_arr);
% para.num_px_thres=length(sim_dist_arr);
para.prominence=prominence;
para.aspect_ratio=aspect_ratio;
para.pixelsize=pixelsize;
para.stepsize=step_size;
para.ADU_thres=ADU_thres;
% kernel argin(mx data, struct(ADU_thres,Hbox_size,Hbox_size_axi,filter_Gau,sim_median_dist,sim_mean_dist,num_px_thres))
%%
if isfield(obj.guihandles.PSFdata,'slicepS')
    obj.guihandles.PSFdata.data=cell(length(path_ary)*obj.guihandles.PSFdata.totslices/obj.guihandles.PSFdata.slicepS,1);
end
for i=1:length(path_ary)
%     [pageNumber ~]=size(imfinfo(path_ary{i,1}));
    imgobj = Tiff(path_ary{i,1},'r');
    if isfield(obj.guihandles.PSFdata,'slicepS')
        tic;
        stack_num=obj.guihandles.PSFdata.totslices/obj.guihandles.PSFdata.slicepS;
        slicepS=obj.guihandles.PSFdata.slicepS;
        data_temp=cell(stack_num,1);
%         for j=(1+slicepS*5):(stack_num+5)*slicepS
%             imgobj.setDirectory(j)
%             temp=imgobj.read;
%             img_full(:,:,j-slicepS*5)=temp((ROI_sele.xmin):(ROI_sele.xmax),(ROI_sele.ymin):(ROI_sele.ymax));
%         end
        img_full = tiffreadVolume(path_ary{i,1});
        img_full = img_full(:,:,1+slicepS*5:end);
        img_seg=cell(stack_num,1);
        for j=1:stack_num
            img_seg{j,1}=img_full(:,:,(j-1)*slicepS+1:j*slicepS);
        end
        
        %cur_thres=a*exp(-ib)
        exp_b=log(ADU_thres/(ADU_thres-ADU_thres_diff))/(length(path_ary)-1);
        exp_a=ADU_thres/exp(-exp_b);
        para.ADU_thres=exp_a*exp(-exp_b*i);
        para.Hbox_size_search=Hbox_size_xy-round(Hbox_size_search_diff*(i-1)/(length(path_ary)-1));
        parfor j=1:stack_num
            img3D=double(img_seg{j,1});
            data_temp{j,1}=PSFseg_kernel(img3D,para);
        end
        for j=1:stack_num
            obj.guihandles.PSFdata.data{(i-1)*obj.guihandles.PSFdata.totslices/obj.guihandles.PSFdata.slicepS+j,1}=data_temp{j,1};
        end
        toc;
        disp(strcat(num2str(i),' hyperstacks out of ',num2str(length(path_ary)),' hyperstacks finished'));
    else
        for k=(ROI_sele.zmin):(ROI_sele.zmax)
            imgobj.setDirectory(k)
            temp=imgobj.read;
            img3D(:,:,(k-ROI_sele.zmin+1))=(temp((ROI_sele.xmin):(ROI_sele.xmax),(ROI_sele.ymin):(ROI_sele.ymax)));
        end
        img3D=double(img3D);
        PSFdata=PSFseg_kernel(img3D,para);
        obj.guihandles.PSFdata.data{i,1}=PSFdata;
    end
end
obj.guihandles.menu_analyze_FM_track.UserData.SegPara = obj.guihandles.PSFdata.SegPara;
obj.guihandles.menu_analyze_FM_track.UserData.AlgorPara = para;
end