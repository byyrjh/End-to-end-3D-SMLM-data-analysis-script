function Fcn_FM_track(hObj,event,obj)
% run PSF extract first on 2 hyperstack trial data 
if isfield(obj.guihandles.PSFdata,'slicepS')
    para = obj.guihandles.menu_analyze_FM_track.UserData.AlgorPara;
    ROI_sele=obj.guihandles.menu_analyze_ROIselect_flag.old;
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
%     path_ary=obj.guihandles.PSFdata.path;
    obj.guihandles.PSFdata.data=cell(length(path_ary)*obj.guihandles.PSFdata.totslices/obj.guihandles.PSFdata.slicepS,1);
    imgobj = Tiff(path_ary{1,1},'r');
    stack_num=obj.guihandles.PSFdata.totslices/obj.guihandles.PSFdata.slicepS;
    slicepS=obj.guihandles.PSFdata.slicepS;
    obj.guihandles.PSFdata.slices=ROI_sele.zmax-ROI_sele.zmin+1;
    obj.guihandles.PSFdata.SegPara = obj.guihandles.menu_analyze_FM_track.UserData.SegPara;
    for j=slicepS*5+1:slicepS*5+slicepS % kickout first 5 stacks
        imgobj.setDirectory(j)
        temp=imgobj.read;
        seed_stack(:,:,j-slicepS*5)=temp((ROI_sele.xmin):(ROI_sele.xmax),(ROI_sele.ymin):(ROI_sele.ymax));
    end
    seed_stack = double(seed_stack);
    seed_seg=PSFseg_kernel(seed_stack,para);
    last_seg.PosSet = seed_seg.PosSet;
    last_seg.VolumeSet = seed_seg.VolumeSet;
    for i=1:length(path_ary)
        tic;
        imgobj = Tiff(path_ary{i,1},'r');
        for j=(1+slicepS*5):(stack_num+5)*slicepS
            imgobj.setDirectory(j)
            temp=imgobj.read;
            img_full(:,:,j-slicepS*5)=temp((ROI_sele.xmin):(ROI_sele.xmax),(ROI_sele.ymin):(ROI_sele.ymax));
        end
        data_raw = cell(stack_num,1);
        data_seg = cell(stack_num,1);
        for j = 1:stack_num
            data_raw{j,1} = double(img_full(:,:,(j-1)*slicepS+1:j*slicepS));
        end
        parfor j=1:stack_num
            data_seg{j,1} = FM_track_kernel(data_raw{j,1},para,last_seg);
        end
        for j = 1:stack_num
            obj.guihandles.PSFdata.data{(i-1)*obj.guihandles.PSFdata.totslices/obj.guihandles.PSFdata.slicepS+j,1}=data_seg{j,1};
        end
        last_seg.PosSet = data_seg{1,1}.PosSet;
        last_seg.VolumeSet = data_seg{1,1}.VolumeSet;
        toc;
        disp(strcat(num2str(i),' hyperstacks out of ',num2str(length(path_ary)),' hyperstacks finished'));
    end
else
    msgbox('This feature is not applicable!');
end
end