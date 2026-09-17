function xyslider_PSFviewer(ObjH,EventData,obj)
hSlider=EventData.AffectedObject;
hPanel=get(hSlider,'parent');
panelname=PanelSelection(obj,hPanel);
index=round(get(panelname.xyslider,'Value'));
panelname.curslice.String=string(index);
axes(panelname.xyshow);
PSFdata=obj.guihandles.PSFdata.data{index,1};
image2disp=PSFdata.xyMIP;
[ims_x ims_y]=size(image2disp);
text_ofs_x=ims_x/100+5;
text_ofs_y=ims_y/100+5;
%% manually define constrast
imagesc(image2disp,[0,100]);
% imagesc(image2disp);
%%
colormap('hot');
axis equal
axis([0 size(image2disp,2) 0 size(image2disp,1)])
hold on
psf_pos=PSFdata.PosSet;
scatter(psf_pos(:,2),psf_pos(:,1),15,'MarkerEdgeColor',[0 1 0]);
if length(panelname.xyshow.Children)>2
    delete(panelname.xyshow.Children(3:4));
end
%%
handleList=allchild(obj.Panel_4.panel);
delete(handleList);
obj.Panel_4.xyshow=axes('Parent',obj.Panel_4.panel);
set(obj.Panel_4.xyshow,'Visible','off');
obj.Panel_4.xyslider=uicontrol('Parent',obj.Panel_4.panel,...
    'Units','normalized','Position',[0.1,0.01,0.8,0.05],...
    'Style','slider',...
    'UserData',struct('row',0,'col',0,'slice',0,'threeDdata',[],'segData',[]));
slider_range=(obj.guihandles.PSFdata.slices);
set(obj.Panel_4.xyslider,'Min',1,'Max',slider_range,'Value',1,'SliderStep',[1/(slider_range-1),1/(slider_range-1)*10]);
obj.Panel_4.xyslider.UserData.slice=slider_range;
obj.Panel_4.curslice=uicontrol(obj.Panel_4.panel,'Style','edit','Units',...
    'normalized','String','1','FontSize',12,'Position',[0.01,0.015,0.08,0.05],...
    'Callback',@(ObjH,EventData) gotoslice(ObjH,EventData,obj));
obj.Panel_4.totslice=uicontrol(obj.Panel_4.panel,'Style','Text','Units',...
    'normalized','String',num2str(slider_range),'FontSize',12,'Position',[0.91,0.015,0.08,0.05]);
obj.Panel_4.xyslider.addlistener('Value','PostSet',@(ObjH,EventData) xyslider_PSFviewer_3D(ObjH,EventData,obj));
%% create position annotation and table
data_temp=cell(size(psf_pos,1),6);
[x y z]=size(PSFdata.VolumeSet{1,1});
center_pos_xy=round(x/2);
center_pos_z=round(z/2);
if isfield(obj.guihandles.PSFdata,'ChiSq')
    idx = obj.guihandles.PSFdata.map_t == index;
    ChiSq = round(obj.guihandles.PSFdata.ChiSq(idx));
    crlb = round(sqrt(obj.guihandles.PSFdata.crlb_res(6,idx))*10); % estimation of localization uncertainty in nm
    LS_pos = round(obj.guihandles.PSFdata.fitting_res(6,idx)*10); % 10 nm stepsize of calibrated PSF in z
    N_photon = round(obj.guihandles.PSFdata.fitting_res(4,idx)*30); % Num of photons conversion factor
    for i=1:size(psf_pos,1)
        data_temp{i,1}=LS_pos(i);
        data_temp{i,2}=crlb(i);
        data_temp{i,3}=ChiSq(i);
        data_temp{i,4}=N_photon(i);
        data_temp{i,5}=PSFdata.PosSet(i,3);
        data_temp{i,6}=any(PSFdata.Manually_sele_idx==i);
    end
else
    for i=1:size(psf_pos,1)
        % text(psf_pos(i,2)-text_ofs_x,psf_pos(i,1)-text_ofs_y,num2str(i),'Color','white','FontSize',10);
        data_temp{i,1}=PSFdata.FWHMset(i,1);
        data_temp{i,2}=PSFdata.FWHMset(i,2);
        data_temp{i,3}=PSFdata.FWHMset(i,3);
        temp=PSFdata.VolumeSet{i,1};
        data_temp{i,4}=temp(center_pos_xy,center_pos_xy,center_pos_z);
        data_temp{i,5}=PSFdata.PosSet(i,3);
        data_temp{i,6}=any(PSFdata.Manually_sele_idx==i);
    end

end
obj.Panel_3.FWHMdisp=uitable(obj.Panel_3.panel,'ColumnName',{'Pos_LS','CRLB_LS','ChiSq','Nphotons',...
    strcat('pos_z_',num2str(obj.guihandles.PSFdata.slices),'/',num2str(center_pos_z)),'rule out'},...
    'Units','normalized','Position',[0 0.3 1 0.6],'FontSize',10);
obj.Panel_3.FWHMdisp.Data=data_temp;
obj.Panel_3.FWHMdisp.ColumnEditable=true;
obj.Panel_3.FWHMdisp.ColumnWidth={55,55,55,60,75,48};
obj.Panel_3.manulsel=uicontrol(obj.Panel_3.panel,'style','pushbutton','Units','normalized','Position',...
    [0.35 0.1 0.3 0.1],'FontSize',10,'String','confirm elimination','Callback',@(hObj,event) PSFruleout(hObj,event,obj,index));
end