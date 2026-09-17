function PSFdata=PSFseg_kernel(img_raw,para)
Gx_oprt=[-1 0 1];
Gy_oprt=[-1;0;1];
Gz_oprt=[-1;0;1];
filter_Gau=para.filter_Gau;
img_median=zeros(size(img_raw));
for i=1:size(img_raw,3)
    img_median(:,:,i)=medfilt2(img_raw(:,:,i));
end
img_filter=imfilter(img_median,filter_Gau);
img_filter_div_=double(div_img(img_filter,Gx_oprt,Gy_oprt,Gz_oprt));
%%
img_filter_div_=imfilter(img_filter_div_,filter_Gau);
%%
%% find local maximum
[x_raw y_raw]=size(img_filter_div_(:,:,1));
img3D_uncrop=double(img_filter_div_(5:x_raw-4,5:y_raw-4,:));
emitter=[];
% ADU_pre_Hmax=0;
Hbox_size_axi=para.Hbox_size_axi;
height_seg=para.height_seg;
width_seg=para.Hbox_size;
Hbox_size=para.Hbox_size_seg; % hard coded as 21*21 pixels 
[x y z]=size(img3D_uncrop);
ADU_thres=para.ADU_thres;
prominence=para.prominence;  % normally prominence is designed to be 0.5
% ADU_cur_Hmax=max(img3D_uncrop(:))*prominence;
% aspect_ratio=para.aspect_ratio;
pixelsize=para.pixelsize;
stepsize=para.stepsize;
[x_,y_,z_]=size(img3D_uncrop);
xx=linspace(1,y_,y_);
yy=linspace(1,x_,x_);
zz=linspace(1,z_,z_);
[xx yy zz]=meshgrid(xx,yy,zz);
% while ADU_cur_Hmax>ADU_thres*prominence
    % setting search crossection
    ADU_cur_Hmax=ADU_thres*prominence;
    idx_lin_uncrop=find(img3D_uncrop>ADU_cur_Hmax);
    peak_value_uncrop=img3D_uncrop(idx_lin_uncrop);
    [row_uncrop,col_uncrop,page_uncrop]=ind2sub(size(img3D_uncrop),idx_lin_uncrop);
    crop_vec=row_uncrop<Hbox_size+1|row_uncrop>x-Hbox_size|col_uncrop<Hbox_size+1|col_uncrop>y-Hbox_size...
        |page_uncrop<Hbox_size_axi+1|page_uncrop>z-Hbox_size_axi;
    peak_value_crop=double(peak_value_uncrop).*double(~crop_vec);
    pixel_map_uncrop=zeros(size(img3D_uncrop)); 
    pixel_map_uncrop(idx_lin_uncrop)=1; 
    idx_map=zeros(size(img3D_uncrop));
    for k=1:length(peak_value_crop)
        idx_map(row_uncrop(k),col_uncrop(k),page_uncrop(k))=k;
    end
    emitter_temp=[];
    while sum(peak_value_crop)>0
        % accept or reject PSF on search hyperplane
        [cur_peak pos]=max(peak_value_crop);
        row_temp=row_uncrop(pos);   % coordinate in uncrop image
        col_temp=col_uncrop(pos);   % cur_sub_vol is ensured to be in cropped image
        page_temp=page_uncrop(pos);
        box_local_max=img3D_uncrop(row_temp,col_temp,page_temp);
        box_candi=img3D_uncrop(row_temp-Hbox_size:row_temp+Hbox_size,col_temp-Hbox_size:col_temp+Hbox_size,page_temp-Hbox_size_axi:page_temp+Hbox_size_axi);
        idx_above_fwhm=find(box_candi>box_local_max*prominence);
        cur_sub_vol=pixel_map_uncrop(row_temp-Hbox_size:row_temp+Hbox_size,col_temp-Hbox_size:col_temp+Hbox_size,page_temp-Hbox_size_axi:page_temp+Hbox_size_axi);
        cur_sub_idx_vol=idx_map(row_temp-Hbox_size:row_temp+Hbox_size,col_temp-Hbox_size:col_temp+Hbox_size,page_temp-Hbox_size_axi:page_temp+Hbox_size_axi);
        lin_idx_sub=find(cur_sub_vol==1);
        peak_value_crop(cur_sub_idx_vol(lin_idx_sub))=0; 
        num_vec_dist=length(idx_above_fwhm)-1;
        if (num_vec_dist>=para.num_px_thres)%&&vec_dist_std<=2*para.sim_dist_std
            %img3D_uncrop(row_temp-Hbox_size:row_temp+Hbox_size,col_temp-Hbox_size:col_temp+Hbox_size,page_temp-Hbox_size_axi:page_temp+Hbox_size_axi)=0;
            sigma_xy=1/sqrt(-2*log(1-1/img3D_uncrop(row_temp,col_temp,page_temp)));
            sigma_z=2/sqrt(-2*log(1-1/img3D_uncrop(row_temp,col_temp,page_temp)));
            damp_factor=(1-exp(-(xx-col_temp).^2/2/sigma_xy^2).*exp(-(yy-row_temp).^2/2/sigma_xy^2).*exp(-(zz-page_temp).^2/2/sigma_z^2));
            img3D_uncrop=img3D_uncrop.*damp_factor;
            emitter_temp=[emitter_temp;row_temp col_temp page_temp];
%         elseif cur_peak*prominence>=0.9*ADU_pre_Hmax
%             img3D_uncrop(row_temp:row_temp,col_temp:col_temp,page_temp:page_temp)=median(reshape(img3D_uncrop(row_temp-1:row_temp+1,col_temp-1:col_temp+1,page_temp-1:page_temp+1),[27 1]));        
        end
    end
    emitter=[emitter;emitter_temp];
%     ADU_pre_Hmax=ADU_cur_Hmax;
%     img_crop=img3D_uncrop(Hbox_size+1:x-Hbox_size,Hbox_size+1:y-Hbox_size,Hbox_size_axi+1:z-Hbox_size_axi);
%     ADU_cur_Hmax=max(img_crop(:))*prominence;
% end
%% fine localization
% noise pixels near real emitter might be identified
% use median filter smooth out noisy pixels
img_bg_sub_fine=uint16(img_median(5:x_raw-4,5:y_raw-4,:));
img_bg_sub_fine=img_bg_sub_fine-median(img_bg_sub_fine,'all');
emitter_fine_row=zeros(size(emitter,1),1);
emitter_fine_col=zeros(size(emitter,1),1);
emitter_fine_page=zeros(size(emitter,1),1);
for p=1:size(emitter,1)
    sub_volume=img_bg_sub_fine(emitter(p,1)-Hbox_size:emitter(p,1)+...
        Hbox_size,emitter(p,2)-Hbox_size:emitter(p,2)+Hbox_size,emitter(p,3)-Hbox_size_axi:emitter(p,3)+Hbox_size_axi);
    idx=find(sub_volume==max(sub_volume(:)));
    [row,col,page]=ind2sub(size(sub_volume),idx);
    emitter_fine_row(p)=emitter(p,1)+row(1)-Hbox_size-1+4;  % +4 converts coordinate to raw untrimed data
    emitter_fine_col(p)=emitter(p,2)+col(1)-Hbox_size-1+4;
    emitter_fine_page(p)=emitter(p,3)+page(1)-Hbox_size_axi-1;
end
[x y z]=size(img_median);
emitter_fine=[emitter_fine_row emitter_fine_col emitter_fine_page];
emitter_fine(emitter_fine(:,3)<height_seg+1,:)=[];
emitter_fine(emitter_fine(:,3)>z-height_seg,:)=[];
emitter_fine(emitter_fine(:,1)<width_seg+1,:)=[];
emitter_fine(emitter_fine(:,1)>x-width_seg,:)=[];
emitter_fine(emitter_fine(:,2)<width_seg+1,:)=[];
emitter_fine(emitter_fine(:,2)>y-width_seg,:)=[];
% emitter_fine=emitter;
%% calculate FWHM
xylength=(width_seg*2+1)*pixelsize;
zlength=(height_seg*2+1)*stepsize;
intp_xy=width_seg*16;
intp_z=height_seg*16;
numofbeads = size(emitter_fine,1);
PSFdata=struct('VolumeSet',[],'xyMIP',[],'PosSet',[],'FWHMset',[],'Manually_sele_idx',[]);
volumeset=cell(numofbeads,1);
measurament_beads=zeros(numofbeads,3);
img_bg_sub_fine=double(img_raw);
for n = 1:numofbeads
    bead = (img_bg_sub_fine(emitter_fine(n,1)-width_seg:emitter_fine(n,1)+width_seg,...
        emitter_fine(n,2)-width_seg:emitter_fine(n,2)+width_seg,...
        emitter_fine(n,3)-height_seg:emitter_fine(n,3)+height_seg));
    volumeset{n,1}=uint16(bead);
    bead=bead/bead(width_seg+1,width_seg+1,height_seg+1);
    try
    measurament_beads(n,1)=FWHM(bead(width_seg+1,:,height_seg+1),xylength,intp_xy);
    measurament_beads(n,2)=FWHM(bead(:,width_seg+1,height_seg+1),xylength,intp_xy);
    measurament_beads(n,3)=FWHM(reshape(bead(width_seg+1,width_seg+1,:),[1,(height_seg*2+1)]),zlength,intp_z);
    catch
        continue
    end
end
measurament_beads(measurament_beads(:,3)==0,:)=[];
measurament_beads(measurament_beads(:,2)==0,:)=[];
measurament_beads(measurament_beads(:,1)==0,:)=[];
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
IMG_MIP=uint16(img_median-median(img_median,'all'));
[~,~,slice_num]=size(IMG_MIP);
IMG_MIP=IMG_MIP(:,:,height_seg:slice_num-height_seg+1);
PSFdata.xyMIP=max(IMG_MIP,[],3);
PSFdata.PosSet=emitter_fine;
PSFdata.FWHMset=measurament_beads;
end