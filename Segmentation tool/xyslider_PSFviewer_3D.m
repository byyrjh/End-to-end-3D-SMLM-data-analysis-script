function xyslider_PSFviewer_3D(ObjH,EventData,obj)
hSlider=EventData.AffectedObject;
hPanel=get(hSlider,'parent');
panelname=PanelSelection(obj,hPanel);
index=round(get(panelname.xyslider,'Value'));
panelname.curslice.String=string(index);
if isempty(obj.Panel_4.xyslider.UserData.threeDdata)
    idx_temp=round(obj.Panel_1.xyslider.Value);
    hyper_stack_idx=ceil(idx_temp/(length(obj.guihandles.PSFdata.data)/length(obj.guihandles.PSFdata.path)));
    stack_idx=idx_temp-(hyper_stack_idx-1)*(length(obj.guihandles.PSFdata.data)/length(obj.guihandles.PSFdata.path));
    path=obj.guihandles.PSFdata.path{hyper_stack_idx,1};
    obj.Panel_4.xyslider.UserData.threeDdata=threeDimgread(path,stack_idx+5,hSlider.Max);
    obj.Panel_4.xyslider.UserData.segData=obj.guihandles.PSFdata.data{idx_temp,1};
end
axes(panelname.xyshow);
image2disp=obj.Panel_4.xyslider.UserData.threeDdata(:,:,index);
%% manually define constrast 
delete(panelname.xyshow.Children);
imagesc(image2disp,[100,200]);  
[ims_x ims_y]=size(image2disp);
text_ofs_x=ims_x/100+5;
text_ofs_y=ims_y/100+5;
% imagesc(image2disp);  
%%
colormap('hot');
axis equal
axis([0 size(image2disp,2) 0 size(image2disp,1)])
hold on
psf_pos=obj.Panel_4.xyslider.UserData.segData.PosSet;
num_emitter = size(psf_pos,1);
num = linspace(1,num_emitter,num_emitter);
num = num(psf_pos(:,3)==index);
pos_arr = psf_pos(psf_pos(:,3)==index,:);
scatter(psf_pos(psf_pos(:,3)==index,2),psf_pos(psf_pos(:,3)==index,1),60,'MarkerEdgeColor',[0 1 0]);
for i = 1:length(num)
    text(pos_arr(i,2)-text_ofs_x,pos_arr(i,1)-text_ofs_y,num2str(num(i)),'Color','white','FontSize',10);
end
end