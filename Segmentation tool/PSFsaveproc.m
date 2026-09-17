function PSFsaveproc(hObj,event,obj)
pararead=findobj('Tag','psfsave');
savepath=findobj('Tag','savePath');
savepath=savepath.String;
Alldata=findobj('Tag','Alldata');
Volumedata=findobj('Tag','Volumedata');
%  % % labels old save routine
% % f=waitbar(0,'1','Name','PSF Data Saving',...
% %     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% % setappdata(f,'canceling',0);

if logical(Alldata.Value)
    data=obj.guihandles.PSFdata;
    save(strcat(savepath,'.mat'),'data','-v7.3');
    obj.guihandles.menu_file_f4.UserData=obj.guihandles.PSFdata.SegPara;
end
if logical(Volumedata.Value)
    data_raw=obj.guihandles.PSFdata.data;
    % %     [aver_psf aver_psf_N]=psf_aver(obj.guihandles.PSFdata.data);
    % %     [~,~,slice]=size(aver_psf);
    % %     for i=1:slice
    % %         imwrite(aver_psf(:,:,i),strcat(savepath,'.tiff'),'WriteMode','append');
    % %         completed=strcat(num2str(round(i/slice,2)*100),'%completed');
    % %         waitbar(i/slice,f,completed);
    % %     end
    num_hypstc=size(data_raw,1);
    num_segs=[];
    for i=1:num_hypstc
        num_segs=[num_segs,size(data_raw{i,1}.VolumeSet,1)];
    end
    [seg_x,seg_y,seg_z]=size(data_raw{1,1}.VolumeSet{1,1});
    seg_sta=zeros(seg_x,seg_y,seg_z*sum(num_segs));
    seg_pos=zeros(sum(num_segs),3);
    for i=1:num_hypstc
        data_raw_temp=data_raw{i,1};
        num_vol_temp=size(data_raw_temp.VolumeSet,1);
        init_z=sum(num_segs(1:i-1))*seg_z;
        init_pos=sum(num_segs(1:i-1));
        seg_pos(init_pos+1:init_pos+num_vol_temp,:)=data_raw_temp.PosSet;
        for j=1:num_vol_temp
            seg_sta(:,:,init_z+1+(j-1)*seg_z:init_z+j*seg_z)=data_raw_temp.VolumeSet{j,1};
        end
    end
    seg_sta=uint16(seg_sta);
    option=struct('message',false);
    res = saveastiff(seg_sta, strcat(savepath,'_seg.tif'),option);
    seg_sta=single(seg_sta);
    loc_pos_x=single(seg_pos(:,1));
    loc_pos_y=single(seg_pos(:,2));
    loc_pos_z=single(seg_pos(:,3));
    save(strcat(savepath,'_seg.mat'),'seg_sta');
    save(strcat(savepath,'_loc_pos_x.mat'),'loc_pos_x');
    save(strcat(savepath,'_loc_pos_y.mat'),'loc_pos_y');
    save(strcat(savepath,'_loc_pos_z.mat'),'loc_pos_z');
end
% delete(f);
close(pararead);
end