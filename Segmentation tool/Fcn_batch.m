function Fcn_batch(hObj,event,obj,imgfile,imgpath)
gpudev=gpuDevice;
num_stack=findobj('Tag','num_stack');
t_wait=findobj('Tag','t_wait');
para_window=findobj('Tag','ParaBatch');
t_wait=str2double(t_wait.String);
num_stack=str2double(num_stack.String);
close(para_window);


%%%%%%
ROI_info=struct('two_ch_corr',[],'xmin',[],'xmax',[],'ymin',[],'ymax',[],'zmin',[],'zmax',[],'coordinate',[],'save_zmin',[],'save_zmax',[]);
ROI_info.two_ch_corr=obj.guihandles.twoChannelcor;
ROI_info.xmin=obj.guihandles.menu_analyze_ROIselect_flag.old.xmin;
ROI_info.xmax=obj.guihandles.menu_analyze_ROIselect_flag.old.xmax;
ROI_info.ymin=obj.guihandles.menu_analyze_ROIselect_flag.old.ymin;
ROI_info.ymax=obj.guihandles.menu_analyze_ROIselect_flag.old.ymax;
ROI_info.zmin=obj.guihandles.menu_analyze_ROIselect_flag.old.zmin;
ROI_info.zmax=obj.guihandles.menu_analyze_ROIselect_flag.old.zmax;
ROI_info.save_zmin=obj.guihandles.saveZrange(1);
ROI_info.save_zmax=obj.guihandles.saveZrange(2);
ROI_info.coordinate='Lab';

median_value=zeros(num_stack,length(imgfile));
bg_value=zeros(num_stack,length(imgfile));
noise_value=zeros(num_stack,length(imgfile));
drifting_Mat=[0 0];
drifting_Mat_accum=[0 0 0];
IniStack=[];
LastStack=[];

slash_pos=strfind(imgfile{1,1},'\');
dot_pos=strfind(imgfile{1,1},'.');
temp_n=length(slash_pos);

imagename=cell(length(imgfile),1);
ini_ind=zeros(length(imgfile),1);
for i=1:length(imgfile)
    imgname=imgfile{i,1}((slash_pos(temp_n)+1):(dot_pos-1));
    imgind=strfind(imgname,'_det');
    ini_ind(i,1)=str2double(imgname(imgind-4:imgind-1));
    imagename{i,1}=imgname;
end

for i=1:num_stack
    for m=1:length(imagename)  % m=1 select one image stack, m=2 select two image stacks
        indx=ini_ind(m,1)-1+1*i;  % indx=ini_ind(m,1)-2+2*i;  Single camera at one time
        indx_4=floor(indx/1000);
        indx_3=floor((indx-1000*indx_4)/100);
        indx_2=floor((indx-1000*indx_4-100*indx_3)/10);
        indx_1=indx-1000*indx_4-100*indx_3-10*indx_2;
        indx_str=strcat(num2str(indx_4),num2str(indx_3),num2str(indx_2),num2str(indx_1));
        temp=imagename{m,1};
        cur_imgname=strcat(temp(1:(imgind-5)),indx_str,temp(imgind:(imgind+4)));
        img_path=strcat(imgpath,cur_imgname,'.tif');
        temp=cell(1,1);
        temp{1,1}=img_path;
        waitflag=0;
        breakflag=false;
        while ~isfile(img_path)
            if waitflag>1
                breakflag=true;
                break
            end
            pause(t_wait);
            waitflag=waitflag+1;
        end
        if breakflag
            continue
        end
        [img_raw breakflag]= tiff_reader_rui(temp);
        
        %%%%%%%%%% constant background subtraction new version %%%%%%%%%%%%
        bg_method=obj.guihandles.batchpipeline.bg_method;
        switch bg_method
            case 'hybrid'
                if isempty(strfind(cur_imgname,'det2'))
                    background=obj.guihandles.batchpipeline.det1;
                else
                    background=obj.guihandles.batchpipeline.det2;
                end
                img=img_raw-background;
                
            case 'camera'
                if isempty(strfind(cur_imgname,'det2'))
                    offset=obj.guihandles.batchpipeline.offset_det1;
                    gain=obj.guihandles.batchpipeline.gain_det1;
                else
                    offset=obj.guihandles.batchpipeline.offset_det2;
                    gain=obj.guihandles.batchpipeline.gain_det2;
                end
                img=uint16((double(img_raw)-offset)./gain);
            otherwise
                img=img_raw;
        end
        % image crop
        if isempty(strfind(cur_imgname,'det2'))
            xmin=obj.guihandles.menu_analyze_ROIselect_flag.old.xmin+obj.guihandles.twoChannelcor(1);
            xmax=obj.guihandles.menu_analyze_ROIselect_flag.old.xmax+obj.guihandles.twoChannelcor(1);
            ymin=obj.guihandles.menu_analyze_ROIselect_flag.old.ymin+obj.guihandles.twoChannelcor(2);
            ymax=obj.guihandles.menu_analyze_ROIselect_flag.old.ymax+obj.guihandles.twoChannelcor(2);
        else
            xmin=obj.guihandles.menu_analyze_ROIselect_flag.old.xmin;
            xmax=obj.guihandles.menu_analyze_ROIselect_flag.old.xmax;
            ymin=obj.guihandles.menu_analyze_ROIselect_flag.old.ymin;
            ymax=obj.guihandles.menu_analyze_ROIselect_flag.old.ymax;
        end
        zmin=obj.guihandles.menu_analyze_ROIselect_flag.old.zmin;
        zmax=obj.guihandles.menu_analyze_ROIselect_flag.old.zmax;
        Image3D=img(xmin:xmax,ymin:ymax,zmin:zmax);
        im_w_bg=img_raw(xmin:xmax,ymin:ymax,zmin:zmax);
        
        % direct current background subtracting
        [x y z]=size(im_w_bg);
        greyval_tot=reshape(im_w_bg,[x*y*z 1]);
        %%%%%%%%%%Gaussian fit%%%%%%%%%%%%%
        greylev_tot=linspace(70,3000,2931);
        [counts_tot centers_tot]=histcounts(greyval_tot,greylev_tot);
        bg_aver=centers_tot(find(counts_tot==max(counts_tot)));
        bg_vec=greyval_tot(find(greyval_tot<(bg_aver+1)));
        greylev=linspace(70,bg_aver,bg_aver-70+1);
        [counts centers]=hist(bg_vec,greylev);
        greylev_sym=linspace(70,70+2*(bg_aver-70),2*(bg_aver-70)+1);
        temp=flip(counts);
        greyval_vec_sym=[counts,temp(2:bg_aver-70+1)];
        [sigma,mu,A]=mygaussfit(greylev_sym,greyval_vec_sym);
        %   f(x) =  a1*exp(-((x-b1)/c1)^2)
        %%%%%%%%%%%%%% maximum and median find %%%%%%%%%%%%%%%%
        [counts_tot centers_tot]=histcounts(greyval_tot,greylev_tot);
        counts_tot(1:(2*(bg_aver-70)+1))=counts_tot(1:(2*(bg_aver-70)+1))-greyval_vec_sym;
        counts_tot(find(counts_tot<0))=0;
        centers_tot_shift=linspace(1,2000,2000);
        counts_tot_shift=counts_tot(bg_aver-70+1:bg_aver-70+2000);
        A=ones(length(counts_tot_shift));
        operator_A=tril(A)';
        posibility_greyval=flip(operator_A*flip(counts_tot_shift'));
        signal_median=abs(posibility_greyval-sum(counts_tot_shift)/2);
        signal_median=find(signal_median==min(signal_median));
        
        
        median_value(i,m)=signal_median;
        bg_value(i,m)=bg_aver;
        noise_value(i,m)=sigma*sqrt(2);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(obj.guihandles.menu_file_f1.UserData.bin_method,"Exter.")
            Image3D=Fcn_imagebin(Image3D,obj.guihandles.menu_file_f1.UserData.binsize);
        end
        for j=1:length(obj.guihandles.batchpipeline.method)
            parameter_method=obj.guihandles.batchpipeline.para{1,j};
            switch obj.guihandles.batchpipeline.method{1,j}
                case 'WBC'
                    resolution_inpix=round(0.4/((obj.guihandles.menu_file_f1.UserData.pixelsize)/(obj.guihandles.menu_file_f1.UserData.magnification)*(obj.guihandles.menu_file_f1.UserData.binsize)));
                    cur_imgname=strcat(cur_imgname,'_WBC');
                    Image3D=Fcn_wbc(Image3D,parameter_method,resolution_inpix);
                case 'DECONV'
                    cur_imgname=strcat(cur_imgname,'_DECONV');
                    data_psf_temp=obj.guihandles.menu_file_f2.UserData;
                    if isstruct(data_psf_temp)
                        [data_psf ~]=psf_aver(data_psf_temp.data);
                    else
                        data_psf=data_psf_temp;
                    end
                    [Image3D breakflag]=Fcn_deconv(Image3D,data_psf,parameter_method(1),parameter_method(2),parameter_method(3),obj.guihandles.menu_file_f1.UserData.binsize);
                case 'Bio'
                    cur_imgname=strcat(cur_imgname,'_Bio');
                    imginfo=struct('magnification',[],'pixelsize',[],'stepsize',[]);
                    imginfo.magnification=obj.guihandles.menu_file_f1.UserData.magnification;
                    imginfo.stepsize=obj.guihandles.menu_file_f1.UserData.stepsize;
                    imginfo.pixelsize=obj.guihandles.menu_file_f1.UserData.pixelsize;
                    imginfo.binsize=obj.guihandles.menu_file_f1.UserData.binsize;
                    ROI_info.coordinate='Bio';
                    Image3D=coordinatetrans(Image3D,imginfo);
                case 'SLICE'
                    cur_imgname=strcat(cur_imgname,'.tif');
                    slice_indx=parameter_method{1,1};
                    savepath=parameter_method{1,2};
                    fullname=strcat(savepath,cur_imgname);
                    cur_imgstack=Image3D(:,:,slice_indx(1):slice_indx(2));
                    if ~isempty(LastStack) & m==1    % always correct drifting with respect to first selected image stack
                        reset(gpudev);
                        [dimx dimy dimz]=size(cur_imgstack);
                        cur_imgstack_GPU=gpuArray(cur_imgstack);
                        IniStack_GPU=gpuArray(IniStack);
                        cross_corr_z=ifftn((fftn(IniStack_GPU)).*conj(fftn(cur_imgstack_GPU)));
                        cross_corr_z=gather(cross_corr_z);
                        clear IniStack_GPU;
                        LastStack_GPU=gpuArray(LastStack);
                        cross_corr=ifftn((fftn(LastStack_GPU)).*conj(fftn(cur_imgstack_GPU)));
                        cross_corr=gather(cross_corr);
                        clear cur_imgstack_GPU LastStack_GPU;
                        [M] = max(cross_corr,[],'all');
                        [M_z] = max(cross_corr_z,[],'all');
                        for q=1:dimz
                            if ~isempty(find(cross_corr(:,:,q)==M))
                                [xshift yshift]=find(cross_corr(:,:,q)==M);
                                %                                 zshift=q;
                                break
                            end
                        end
                        
                        for q=1:dimz
                            if ~isempty(find(cross_corr_z(:,:,q)==M_z))
                                zshift=q;
                                break
                            end
                        end
                        
                        if xshift>round(dimx/2)
                            xshift=xshift-dimx;
                        end
                        if yshift>round(dimy/2)
                            yshift=yshift-dimy;
                        end
                        
                        
                        if zshift>round(dimz/2)
                            zshift=zshift-dimz;
                        end
                        
                        zshift=zshift-1;
                        xshift=xshift-1;
                        yshift=yshift-1;
                        drifting_Mat=[drifting_Mat;xshift yshift];
                        drifting_vec_accum=[sum(drifting_Mat,1) zshift];
                        drifting_Mat_accum=[drifting_Mat_accum;drifting_vec_accum];
                    end
                    if m==1
                        LastStack=cur_imgstack;
                    end
                    if isempty(IniStack)
                        IniStack=cur_imgstack;
                    end
                    xshift=drifting_Mat_accum(end,1);
                    yshift=drifting_Mat_accum(end,2);
                    zshift=drifting_Mat_accum(end,3);
                    [dimx dimy]=size(Image3D(:,:,1));
                    datatosave=uint16(zeros(dimx,dimy,slice_indx(2)-slice_indx(1)+1));
                    for k=slice_indx(1)-zshift:slice_indx(2)-zshift
                        datatosave(:,:,k-slice_indx(1)+zshift+1)=uint16(imtranslate(Image3D(:,:,k),[yshift xshift]));
                    end
                    option=struct('message',false);
                    res = saveastiff(datatosave, fullname,option);
            end
        end
        %%%%%%%%%%%%%%% one image stack processing finished %%%%%%%%%%%%%%%
    end
end
SNR=median_value./noise_value;
SBR=median_value./bg_value;
photobleaching=struct('SNR',SNR,'SBR',SBR);
Batch_log=struct('ROI_info',ROI_info,'photobleaching_data',photobleaching,'drifting_Mat_accum',drifting_Mat_accum);
save(strcat(savepath,'Batch_log.mat'),'Batch_log');
end