function Fcn_tifhyperstackread(hObj,event,obj)
obj.guihandles.ifhyperstack=true;
handleList=allchild(obj.Panel_1.panel);
delete(handleList);
handleList=allchild(obj.Panel_2.panel);
delete(handleList);
handleList=allchild(obj.Panel_3.panel);
delete(handleList);
handleList=allchild(obj.Panel_4.panel);
delete(handleList);
obj.guihandles.Data=cell(2,5);
obj.guihandles.batchpipeline.method=[];
obj.guihandles.batchpipeline.para=[];
obj.Panel_1.xyshow=axes('Parent',obj.Panel_1.panel);
set(obj.Panel_1.xyshow,'Visible','off');
obj.Panel_1.xyslider=uicontrol('Parent',obj.Panel_1.panel,...
    'Units','normalized','Position',[0.1,0.01,0.8,0.05],...
    'Style','slider',...
    'UserData',struct('row',0,'col',0,'slice',0,'threeDdata',[],'xyMIP',[]));
obj.Panel_1.curslice=uicontrol(obj.Panel_1.panel,'Style','edit','Units',...
    'normalized','String','0','FontSize',12,'Position',[0.01,0.015,0.08,0.05],...
    'Callback',@(ObjH,EventData) gotoslice(ObjH,EventData,obj));
obj.Panel_1.totslice=uicontrol(obj.Panel_1.panel,'Style','Text','Units',...
    'normalized','String','0','FontSize',12,'Position',[0.91,0.015,0.08,0.05]);
obj.Panel_1.xyslider.addlistener('Value','PostSet',@(ObjH,EventData) xyslider_continousshow(ObjH,EventData,obj));
set(obj.Panel_1.xyslider,'Enable','off');


%%
pararead=figure('Name','image information register','MenuBar','none','ToolBar','none','Tag','Inforeg');
initPosition = pararead.Position;  % initPosition(1) initPosition(2) are x y cordinate, whose origin is located on lower left corner of screen
pararead.Position = [initPosition(1), ...
    initPosition(2) - 400 + initPosition(4), 300, 300];
pixelsize=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','Camera pixel size(micron)','FontSize',10,...
    'Position',[0.05,0.85,0.3,0.1]);
magnification=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','Magnification','FontSize',10,...
    'Position',[0.04,0.65,0.3,0.1]);
step_size=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','Step size(micron)','FontSize',10,...
    'Position',[0.03,0.5,0.4,0.1]);
bin_size=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','Bin size','Fontsize',10,...
    'Position',[0.003,0.35,0.3,0.1]);
ROI_raw=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','ROI_raw','Fontsize',10,...
    'Position',[0.05,0.2,0.2,0.1]);
stepsize=string(obj.guihandles.menu_file_f6.UserData.stepsize);
pixsize=string(obj.guihandles.menu_file_f6.UserData.pixelsize);
magnif=string(obj.guihandles.menu_file_f6.UserData.magnification);
binsize=string(obj.guihandles.menu_file_f6.UserData.binsize);
ROI_L_mem=string(obj.guihandles.menu_file_f6.UserData.roi.L);
ROI_T_mem=string(obj.guihandles.menu_file_f6.UserData.roi.T);
ROI_R_mem=string(obj.guihandles.menu_file_f6.UserData.roi.R);
ROI_B_mem=string(obj.guihandles.menu_file_f6.UserData.roi.B);
pixelsize_in=uicontrol(pararead,'Style','edit','Units','normalized',...
    'Tag','pixelsize','String',pixsize,'FontSize',10,'Position',[0.55,0.84,0.3,0.1]);
magnification_in=uicontrol(pararead,'Style','edit','Units','normalized',...
    'Tag','magnification','String',magnif,'FontSize',10,'Position',[0.55,0.67,0.3,0.1]);
step_size_in=uicontrol(pararead,'Style','edit','Units','normalized',...
    'Tag','step_size','FontSize',10,'String',stepsize,'Position',[0.55,0.52,0.3,0.1]);
bin_size_in=uicontrol(pararead,'Style','edit','Units','normalized',...
    'Tag','bin_size','FontSize',10,'String',binsize,'Position',[0.25,0.37,0.1,0.1]);
ROI_L=uicontrol(pararead,'Style','edit','Units','normalized',...
    'Tag','ROI_L','FontSize',10,'String',ROI_L_mem,'Position',[0.35,0.22,0.1,0.1]);
ROI_T=uicontrol(pararead,'Style','edit','Units','normalized',...
    'Tag','ROI_T','FontSize',10,'String',ROI_T_mem,'Position',[0.52,0.22,0.1,0.1]);
ROI_R=uicontrol(pararead,'Style','edit','Units','normalized',...
    'Tag','ROI_R','FontSize',10,'String',ROI_R_mem,'Position',[0.69,0.22,0.1,0.1]);
ROI_B=uicontrol(pararead,'Style','edit','Units','normalized',...
    'Tag','ROI_B','FontSize',10,'String',ROI_B_mem,'Position',[0.86,0.22,0.1,0.1]);
ROI_L_label=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','L','FontSize',10,...
    'Position',[0.3,0.25,0.05,0.05]);
ROI_T_label=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','T','FontSize',10,...
    'Position',[0.46,0.25,0.05,0.05]);
ROI_R_label=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','R','FontSize',10,...
    'Position',[0.63,0.25,0.05,0.05]);
ROI_B_label=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','B','FontSize',10,...
    'Position',[0.8,0.25,0.05,0.05]);


binmethod=uicontrol(pararead,'Style','text','Units','normalized',...
    'String','Method','FontSize',10,...
    'Position',[0.35,0.35,0.2,0.1]);
binsele_in=uibuttongroup(pararead,'Visible','on','Position',[0.53,0.37,0.45,0.1],'Tag','bin_sele');
bin_inter=uicontrol(binsele_in,'Style','radiobutton','String','Inter.','Position',[10,1,60,30],...
    'HandleVisibility','on','FontSize',10);
bin_exter=uicontrol(binsele_in,'Style','radiobutton','String','Exter.','Position',[60,1,60,30],...
    'HandleVisibility','on','FontSize',10);

confirmation=uicontrol(pararead,'Style','pushbutton','String','confirm',...
    'Units','normalized','Position',[0.4 0.06 0.2 0.1],'FontSize',10,...
    'Callback',{@obj.Inforegister,obj});

%%
[imgfile, imgpath] = uigetfile_rui( {'*.tiff';'*.tif'}, 'please select image');
obj.guihandles.ifImg=(length(imgfile)==1);
if ~obj.guihandles.ifImg    % these callback function can either read image raw be processed or PSF raw data
    obj.guihandles.PSFdata.path=imgfile;
    obj.guihandles.PSFdata.data=cell(size(imgfile)); % every cell element should be a structure
    temp=cell(1);
    temp{1,1}=imgfile{1,1};
else
    temp=imgfile;
end
%%




[img breakflag]= tiff_reader_rui(temp);
if ~breakflag
    [x y z]=size(img);
    data=struct('row',x,'col',y,'slice',z,'threeDdata',img,'xyMIP',[]);
    slash_pos=strfind(imgfile{1,1},'\');
    dot_pos=strfind(imgfile{1,1},'.');
    temp_n=length(slash_pos);
    imagename=imgfile{1,1};
    imagename=imagename((slash_pos(temp_n)+1):(dot_pos-1));
    obj.guihandles.menu_file_f3.UserData.imagename=imagename;
    if obj.guihandles.ifImg
        set(obj.guihandles.menu_analyze,'Enable','on');
        set(obj.guihandles.handlearr.img,'Enable','on');
        set(obj.guihandles.handlearr.psf,'Enable','off');
        obj.guihandles.menu_file_f6.UserData.row=x;   % original data
        obj.guihandles.menu_file_f6.UserData.col=y;
        obj.guihandles.menu_file_f6.UserData.slice=z;
        obj.guihandles.menu_file_f6.UserData.threeDdata=img;
        obj.guihandles.menu_file_f3.UserData.imagepath=imgpath;
    else
        set(obj.guihandles.menu_analyze,'Enable','on');
        set(obj.guihandles.handlearr.img,'Enable','off');
        set(obj.guihandles.handlearr.psf,'Enable','on');
        obj.guihandles.PSFdata.name=imagename;
    end
    set(obj.Panel_1.xyslider,'Enable','on','Min',1,'Max',z,'Value',1,'SliderStep',[1/z,1/z*10]);
    set(obj.Panel_1.xyslider,'UserData',data);   % data to be analyzed
    obj.Panel_1.totslice.String=string(z);
    axes(obj.Panel_1.xyshow);
    imshow(img(:,:,1),[]);
end

end