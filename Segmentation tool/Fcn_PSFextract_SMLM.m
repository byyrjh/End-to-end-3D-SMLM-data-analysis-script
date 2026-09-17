function Fcn_PSFextract(hObj,event,obj)
pararead=findobj('Tag','PSFpara');
width_in=findobj('Tag','width');
height_in=findobj('Tag','height');
% parameters setting
width=str2double(width_in.String);  % usually 30 with lattice light sheet microscopy configuration
height=str2double(height_in.String);  % unit in step pixel
stepsize=obj.guihandles.PSFdata.stepsize;  % unit um
pixelsize=(obj.guihandles.PSFdata.pixelsize)/(obj.guihandles.PSFdata.magnification);  % unit um
path_ary=obj.guihandles.PSFdata.path;
ROI_sele=obj.guihandles.menu_analyze_ROIselect_flag.old;
obj.guihandles.PSFdata.slices=ROI_sele.zmax-ROI_sele.zmin+1;
close(pararead);

% f=waitbar(0,'0','Name','PSF processing',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% setappdata(f,'canceling',0);
% pages_tot=(ROI_sele.zmax-ROI_sele.zmin+1)*length(path_ary);

%%   new algorithm parameter settings
pixel_size=pixelsize*1000;  % unit nm
step_size=stepsize*1000; % unit nm
filter_size=[5 5 5];
std=1; % 3D Gaussian filter
Gx_oprt=[-1 0 1];
Gy_oprt=[-1;0;1];
Gz_oprt=[-1;0;1];
Hbox_size=width;
Hbox_size_axi=height; % unit in step pixel
num_px_thres=19;  % NA 1.1    lambda 515 nm
ADU_thres=100;
global img3D_proc;
filter_size_in=(filter_size-1)/2;
[xx,yy,zz] = meshgrid(-filter_size_in(2):filter_size_in(2),-filter_size_in(1):filter_size_in(1),-filter_size_in(3):filter_size_in(3));
aspect_ratio=step_size/pixel_size;
filter_Gau=exp(-(xx.^2+yy.^2+(aspect_ratio*zz).^2)/(2*std^2));
filter_Gau(filter_Gau<eps*max(filter_Gau(:))) = 0;
filter_Gau=filter_Gau/sum(filter_Gau(:));

sigma=sqrt(2/log(2));  %FWHM=4 pixels such that exp(-2^2/2/sigma^2)=1/2
sim_x=linspace(-Hbox_size,Hbox_size,Hbox_size*2+1);
[sim_xx sim_yy sim_zz]=meshgrid(sim_x,sim_x,sim_x);
sim_ADU=exp(-(sim_xx.^2+sim_yy.^2+sim_zz.^2)/2/sigma^2);
sim_rr=sqrt(sim_xx.^2+sim_yy.^2+sim_zz.^2);
sim_idx=find(sim_ADU>0.5);
sim_dist_arr=sim_rr(sim_idx);
sim_median_dist=median(sim_dist_arr);  % median is independent of Hbox_size_axi and step size as simulated PSF is isotropic

%%
for i=1:length(path_ary)
    [pageNumber ~]=size(imfinfo(path_ary{i,1}));
    imgobj = Tiff(path_ary{i,1},'r');
    img3D=zeros(ROI_sele.xmax-ROI_sele.xmin+1,ROI_sele.ymax-ROI_sele.ymin+1,ROI_sele.zmax-ROI_sele.zmin+1);
    for j=(ROI_sele.zmin):(ROI_sele.zmax)
        imgobj.setDirectory(j)
        temp=imgobj.read;
        img3D(:,:,(j-ROI_sele.zmin+1))=temp((ROI_sele.xmin):(ROI_sele.xmax),(ROI_sele.ymin):(ROI_sele.ymax));
%         if getappdata(f,'canceling')
%             delete(f);
%             return
%         end
%         page_idx=(i-1)*(ROI_sele.zmax-ROI_sele.zmin+1)+(j+1-ROI_sele.zmin);
%         completed=strcat(num2str(round(page_idx/pages_tot,2)*100),'%completed');
%         waitbar(page_idx/pages_tot,f,completed);
    end
    %%  new algorithm calc
    img_raw=double(img3D);
    % filter raw data
    img_filter=imfilter(img_raw,filter_Gau);
    % calculate Divergence
    img_filter_div=double(div_img(img_filter,Gx_oprt,Gy_oprt,Gz_oprt));
    % filter divergence data
    img_filter_div_filter=imfilter(img_filter_div,filter_Gau);
    % remove edge
    [ims_x ims_y ims_z]=size(img_filter_div_filter);
    img3D_edge_free=uint16(img_filter_div_filter(5:ims_x-4,5:ims_y-4,:));
    
    ADU_cur_Hmax=max(img3D_edge_free(:))/2;
    % protect edge removed raw data
    img3D_proc=img3D_edge_free;
    emitter=[];
    ADU_pre_Hmax=0;
    [x y z]=size(img3D_proc);
    img_crop=img3D_proc(Hbox_size+1:x-Hbox_size,Hbox_size+1:y-Hbox_size,Hbox_size_axi+1:z-Hbox_size_axi);
    %% core part
    while ADU_cur_Hmax>ADU_thres/2
        idx_lin=find(img_crop>ADU_cur_Hmax);
        [row, col, page] = ind2sub(size(img_crop), idx_lin);
        peak_value=img_crop(idx_lin);
        emitter=[emitter;
            loc_core_3D(row,col,page,peak_value,aspect_ratio,num_px_thres,ADU_pre_Hmax,Hbox_size,Hbox_size_axi,sim_median_dist)];
        ADU_pre_Hmax=ADU_cur_Hmax;
        img_crop=img3D_proc(Hbox_size+1:x-Hbox_size,Hbox_size+1:y-Hbox_size,Hbox_size_axi+1:z-Hbox_size_axi);
        ADU_cur_Hmax=max(img_crop(:))/2;
    end
    %% fine localization
    img_bg_sub_fine=uint16(img_raw(5:ims_x-4,5:ims_y-4,:));
    emitter_fine_row=zeros(size(emitter,1),1);
    emitter_fine_col=zeros(size(emitter,1),1);
    emitter_fine_page=zeros(size(emitter,1),1);
    for p=1:size(emitter,1)
        sub_volume=img_bg_sub_fine(emitter(p,1)-Hbox_size:emitter(p,1)+...
            Hbox_size,emitter(p,2)-Hbox_size:emitter(p,2)+Hbox_size,emitter(p,3)-Hbox_size_axi:emitter(p,3)+Hbox_size_axi);
        idx=find(sub_volume==max(sub_volume(:)));
        [row,col,page]=ind2sub(size(sub_volume),idx);
        emitter_fine_row(p)=emitter(p,1)+row(1)-Hbox_size-1;
        emitter_fine_col(p)=emitter(p,2)+col(1)-Hbox_size-1;
        emitter_fine_page(p)=emitter(p,3)+page(1)-Hbox_size_axi-1;
    end
    emitter_fine=[emitter_fine_row emitter_fine_col emitter_fine_page];
    emitter_fine(emitter_fine(:,3)<height+1,:)=[];
    emitter_fine(emitter_fine(:,3)>z-height,:)=[];
    emitter_fine(emitter_fine(:,1)<width+1,:)=[];
    emitter_fine(emitter_fine(:,1)>x-width,:)=[];
    emitter_fine(emitter_fine(:,2)<width+1,:)=[];
    emitter_fine(emitter_fine(:,2)>y-width,:)=[];
    %%
    xylength=(width*2+1)*pixelsize;
    zlength=(height*2+1)*stepsize;
    intp_xy=width*16;
    intp_z=height*16;
    numofbeads = size(emitter_fine,1);
    PSFdata=struct('VolumeSet',[],'xyMIP',[],'PosSet',[],'FWHMset',[],'Manually_sele_idx',[]);
    volumeset=cell(numofbeads,1);
    measurament_beads=zeros(numofbeads,3);
    img_bg_sub_fine=double(img_bg_sub_fine);
    img_bg_sub_fine=img_bg_sub_fine-107;
    for n = 1:numofbeads
        bead = (img_bg_sub_fine(emitter_fine(n,1)-width:emitter_fine(n,1)+width,...
            emitter_fine(n,2)-width:emitter_fine(n,2)+width,...
            emitter_fine(n,3)-height:emitter_fine(n,3)+height));
        volumeset{n,1}=uint16(bead);
        bead=bead/bead(width+1,width+1,height+1);
        measurament_beads(n,1)=FWHM(bead(width+1,:,height+1),xylength,intp_xy);
        measurament_beads(n,2)=FWHM(bead(:,width+1,height+1),xylength,intp_xy);
        measurament_beads(n,3)=FWHM(reshape(bead(width+1,width+1,:),[1,(height*2+1)]),zlength,intp_z);
    end
    inx_min=find(measurament_beads(:,3)==min(measurament_beads(:,3)));
    if length(inx_min)>1
        inx_min=inx_min(1);
    end
    inx_max=find(measurament_beads(:,3)==max(measurament_beads(:,3)));
    if length(inx_max)>1
        inx_max=inx_max(1);
    end
        PSFdata.Manually_sele_idx=[];

%     PSFdata.Manually_sele_idx=[inx_max inx_min];
    PSFdata.VolumeSet=volumeset;
    PSFdata.xyMIP=MIP(img_bg_sub_fine,1,1,1,'xy');
    PSFdata.PosSet=emitter_fine;
    PSFdata.FWHMset=measurament_beads;
    obj.guihandles.PSFdata.data{i,1}=PSFdata;
end
% delete(f);
% volumeViewer(average_bead_8);
end