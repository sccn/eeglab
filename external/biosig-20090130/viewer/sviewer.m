function varargout = sviewer(varargin)
% SVIEWER
% Select HELP in the Info menu 
%
% $Id: sviewer.m,v 1.1 2009-01-30 06:04:53 arno Exp $
% Copyright by (C) 2004 Franz Einspieler <znarfi5@hotmail.com> and
%              (C) 2004,2008 Alois Schloegl <a.schloegl@ieee.org>
% University of Technology Graz, Austria
% This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sviewer_OpeningFcn, ...
                   'gui_OutputFcn',  @sviewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

switch nargin,
 case 0,
  if ~isempty(findobj('Tag', 'sviewer'))
    error('Only one copy of SViewer may be run');
  end
 case 1,
  switch varargin{1}
   case 'draw_detection',
       drawdetection;
   case 'drawline'
       Drawline_Callback;
   case 'Loadfile_Detection'
       Loadfile_Detection;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before sviewer is made visible.
function sviewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sviewer (see VARARGIN)

% Choose default command line output for sviewer
handles.output = hObject;
set(gcf,'Color',[0.949,0.949,1]);
set(gcf,'Name','SViewer (C) 2004,2008 V1.52, University of Technology Graz, Austria');
set(findobj('Tag','ChannelsMenu'),'Enable','off');

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = sviewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OptionMenu_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create menu
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close loaded file
function CloseFile_Callback(hObject, eventdata, handles)
% hObject    handle to NewWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
if isempty(Data)
    errordlg('No File is open!', 'Error');
    return;
else
    if isfield(Data,'HDR')
        Data.HDR = sclose(Data.HDR);
    end
    deleteObj(Data,0);
    Data = [];
end
set(findobj('Tag','checkbox_grid'),'Enable','off');
set(findobj('Tag','checkbox_range'),'Enable','off');
set(findobj('Tag','numb_act_channels'),'Enable','off');
set(findobj('Tag','Slider_ShowChannels'),'Enable','off');
set(findobj('Tag','zoom_factor'),'Enable','off');
set(findobj('Tag','pb_zoom_in_all'),'Enable','off');
set(findobj('Tag','pb_zoom_out_all'),'Enable','off');
set(findobj('Tag','pb_Autoscale'),'Enable','off');
set(findobj('Tag','displayed'),'Enable','off');
set(findobj('Tag', 'goto_h'), 'Enable','off');
set(findobj('Tag', 'goto_min'), 'Enable','off'); 
set(findobj('Tag', 'goto_sec'), 'Enable','off'); 
set(findobj('Tag', 'edit_gotosecond'), 'Enable','off'); 
set(findobj('Tag', 'checkbox_buttons'), 'Enable','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open a new window
function PrintWindow_Callback(hObject, eventdata, handles)
% hObject    handle to NewWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
if isempty(Data)
    errordlg('No File is open!', 'Error');
    return;
else
    printdlg(findobj('Tag', 'sviewer'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open a selected data file
function OpenFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile({ ...
            '*.*', 'All Files (*.*)'; ...
            '*.edf', 'EDF-Files (*.edf)'; ...
            '*.bkr', 'BKR-Files (*.bkr)'; ...
            '*.bdf', 'BDF-Files (*.bdf)'; ...
            '*.cnt', 'CNT-Files (*.cnt)'; ...
            '*.rec', 'REC-Files (*.rec)'; ...
            '*.hea', 'HEA-Files (*.hea)'; ...
            '*.au', 'AU-Files (*.au)'; ...
            '*.wav', 'WAV-Files (*.wav)'; ...
            '*.sig', 'SIG-Files (*.sig)'; ...
            '*.eeg', 'EEG-Files (*.eeg)'; ...
            '*.sma', 'SMA-Files (*.sma)'; ...
            '*.gdf', 'GDF-Files (*.gdf)'}, ...
            'Open file');

% [file,path] = uigetfile({ ...
%                         '*.*', 'All Files (*.*)'}, ...
%                         'Open file');

if file == 0
    return;
else
    Data = get(findobj('Tag','sviewer'),'UserData');
    detect_stat = get(findobj('Tag','Startdetection'),'Label');
    if strmatch('Stop/Save Detection',detect_stat)
        setdefault_detection_stop(Data,1);
    end
        setdefault(file,path);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check and setting default values
function setdefault(file,path)

Data = get(findobj('Tag','sviewer'),'UserData');

if isfield(Data,'HDR')
    try Data.HDR = sclose(Data.HDR); end
    deleteObj(Data,0);
    Data=[];
end

global BIOSIG_GLOBAL;
if ~isfield(BIOSIG_GLOBAL,'ISLOADED_EVENTCODES') || ~BIOSIG_GLOBAL.ISLOADED_EVENTCODES,
	%%% define BIOSIG_GLOBAL.EVENT
	sclose(sopen('eventcodes.txt'));
end;
Data.Eventcodes_txt = BIOSIG_GLOBAL.EVENT;
Event_string = sprintf('%s |',Data.Eventcodes_txt.GroupDesc{:});
Event_string = Event_string(1:end-2);
set(findobj('Tag','Event'),'String',Event_string);
set(findobj('Tag','listbox_Event'),'String',Event_string);

%pos_eventdetail = find(Data.Eventcodes_txt.CodeIndex >= Data.Eventcodes_txt.GroupValue(1) ...
%                            & Data.Eventcodes_txt.CodeIndex < Data.Eventcodes_txt.GroupValue(2));
select = 1; s = Data.Eventcodes_txt;
pos_eventdetail = find(bitand(s.CodeIndex,s.GroupMask(select))==s.GroupValue(select));        % decoding of eventtable 
Eventdetail_string = sprintf('%s |',Data.Eventcodes_txt.CodeDesc{pos_eventdetail});
Eventdetail_string = Eventdetail_string(1:end-2);

set(findobj('Tag','Event_detail'),'String',Eventdetail_string);
set(findobj('Tag','listbox_Event_detail'),'String',['all Events |' Eventdetail_string]);
Data.detcolor = load('detcolor.mat');
set(findobj('Tag', 'Slider1'), 'Value',0);

newfile = [path,file];
Data.HDR = sopen(newfile,'r',0,'OVERFLOWDETECTION:OFF');
Data.Total_length_samples = Data.HDR.NRec * max(Data.HDR.SPR);
Data.Total_length_sec = Data.Total_length_samples / max(Data.HDR.SampleRate);
Data.ShowSamples = min(1000,Data.Total_length_samples);

Data.NS = Data.HDR.NS;

exist_PhysMin = isfield(Data.HDR,'PhysMin');
exist_PhysMax = isfield(Data.HDR,'PhysMax');

if exist_PhysMin == 0
    Data.HDR.PhysMin = [];
    size_Min = 0;
else
    size_Min = size(Data.HDR.PhysMin);
end

if size_Min(1) == 0
    Data.HDR.PhysMin = [ones(1,Data.HDR.NS)*(-999)]';
else
    if size_Min(1) < Data.HDR.NS
        Data.HDR.PhysMin = [ones(1,Data.HDR.NS)*-999]';
    end
end

if exist_PhysMax == 0
    Data.HDR.PhysMax = [];
    size_Max = 0;
else
    size_Max = size(Data.HDR.PhysMax);
end

if size_Max(1) == 0
    Data.HDR.PhysMax = [ones(1,Data.HDR.NS)*(999)]';
else
    if size_Max(1) < Data.HDR.NS
        Data.HDR.PhysMax = abs(Data.HDR.PhysMin);
    end
end

pos_Min_false = find(Data.HDR.PhysMin > 0);
pos_Max_false = find(Data.HDR.PhysMax < 0);
Data.HDR.PhysMin(pos_Min_false,:) = Data.HDR.PhysMin(pos_Min_false,:)*(-1);
Data.HDR.PhysMax(pos_Max_false,:) = Data.HDR.PhysMax(pos_Max_false,:)*(-1);

if ~isfield(Data.HDR,'PhysDim')
    Data.HDR.PhysDim = {''};
end

if ~isfield(Data.HDR,'Label')
        Data.HDR.Label = {};
end

if ischar(Data.HDR.Label)
        Data.HDR.Label = cellstr(Data.HDR.Label);
end;
size_label = length(Data.HDR.Label);
Data.NS_max = Data.HDR.NS;
Data.Channel = cell(Data.HDR.NS,2);

if size_label == 1
        Data.HDR.Label = repmat(Data.HDR.Label,Data.HDR.NS,1);
elseif size_label < Data.HDR.NS,
        for i = size_label+1 : Data.HDR.NS,
                Data.HDR.Label{i} = ['Channel ' int2str(i)];
        end;
end;
Data.Channel = Data.HDR.Label(:);
Data.Channel(:,2) = num2cell([1:Data.HDR.NS]');
Data.allChannel = Data.HDR.Label(:);


Data.ShowChannelmax = 4;
set(findobj('Tag','numb_act_channels'),'String',min(size(Data.Channel,1),Data.ShowChannelmax));
set(findobj('Tag','Slider_Channel'),'Visible','on','Value',1);
Data.ChannelConf.Display_min = Data.HDR.PhysMin;
Data.ChannelConf.Display_max = Data.HDR.PhysMax;
Data.ChannelConf.Scale = 2;
Data.Slider.tsec = 0;
Data.File.file = file;
Data.File.path = path;
set(findobj('Tag','sviewer'),'UserData',Data);

drawnew (file,path,Data.ShowSamples,0); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawnew (file,path,showsamples,only_plot)
newfile = [path,file];
Data = get(findobj('Tag','sviewer'),'UserData');
if ~isequal(Data, 0)
    if isfield(Data,'NoS');
        deleteObj(Data,only_plot);
        Data.NS = [];
    end
    Data = get(findobj('Tag','sviewer'),'UserData');
    Data.NoS = showsamples / max(Data.HDR.SampleRate);
    Data.NS = str2num(get(findobj('Tag','numb_act_channels'),'String'));
    tsec = Data.Total_length_sec;
    tmin = floor(tsec/60);
    tsec = rem(tsec,60);
    th = floor(tmin/60);
    tmin = rem(tmin,60);

    set(findobj('Tag','total_length'), 'String', ...
        sprintf('%02d:%02d:%02d', th,tmin,floor(tsec))); % total length
    
    titlename = get(gcf,'Name');
    oldfilepos = strfind(titlename,':');
    if ~isempty(oldfilepos)
        titlename = deblank(titlename(1:oldfilepos-6));
    end
    titlename = [titlename, '             File: ' file];
    set(gcf,'Name',titlename);
    
    set(findobj('Tag','displayed'), 'String', ...
        sprintf('%c',num2str(Data.NoS))); % displayed length
    actchannels = size(Data.Channel);
    set(findobj('Tag','numb_all_channels'), 'String', Data.HDR.NS);
    set(findobj('Tag','sviewer'),'UserData',Data);
    
    Data = get(findobj('Tag','sviewer'),'UserData');
    
    try
        pos_s1 = get(findobj('Tag', 'Slider1'), 'Value');
        Data.Slider.Pos = pos_s1;
        length_s1 = Data.Total_length_samples;
        tsec_s1 = pos_s1 * (length_s1 - Data.ShowSamples) / max(Data.HDR.SampleRate);
        tmin_ts = floor(tsec_s1 / 60);
        tsec1_ts = rem(tsec_s1,60);
        th_ts = floor(tmin_ts / 60);
        tmin_ts = rem(tmin_ts,60); 
        tsec_s1 = round(tsec_s1 * Data.HDR.SampleRate) / Data.HDR.SampleRate;
    catch
        tsec_s1 = 0;
        Data.Slider.Pos = 0;
        th_ts = 0;
        tmin_ts = 0; 
        tsec1_ts = 0;
    end

    [Data.signal,Data.HDR] = sread(Data.HDR,Data.NoS,tsec_s1);  
    numb_samples = Data.Total_length_samples;
    max_length = numb_samples - Data.ShowSamples;
    step = Data.ShowSamples / max_length / 2; 
    try
        set(findobj('Tag','Slider1'), ...
            'Enable','on', ...
            'Units', 'Normalized', ...
            'Max', 1, ...      
            'Min', 0, ...
            'SliderStep',[step step*2]);
    catch
        errordlg('Incorrect file! (Error:S1)', 'Error');
        CloseFile_Callback;
        return;
    end
    if (size(Data.Channel,1) - Data.NS) <= 0
        set(findobj('Tag','Slider_Channel'),'Enable','off');
    else
        step_channel = 1 / (size(Data.Channel,1) - Data.NS);
        set(findobj('Tag','Slider_Channel'), ...
            'Enable','on', ...
            'Units', 'Normalized', ...
            'Max', 1, ...      
            'Min', 0, ...
            'SliderStep',[step_channel step_channel*10]);
    end
     
    step_showchannel = 1 / size(Data.Channel,1);
    if step_showchannel == inf
        set(findobj('Tag','Slider_ShowChannels'),'Enable','on');
    else
        set(findobj('Tag','Slider_ShowChannels'), ...
            'Enable','on', ...
            'Units', 'Normalized', ...
            'Max', 1, ...      
            'Min', 0, ...
            'Value', step_showchannel * Data.NS, ...
            'SliderStep',[step_showchannel step_showchannel*10]);
    end
       
    % activate checkboxes
    set(findobj('Tag','checkbox_grid'),'Enable','on');
    set(findobj('Tag','checkbox_range'),'Enable','on');
    set(findobj('Tag','numb_act_channels'),'Enable','on');
    set(findobj('Tag','Slider_ShowChannels'),'Enable','on');
    set(findobj('Tag','zoom_factor'),'Enable','on');
    set(findobj('Tag','pb_zoom_in_all'),'Enable','on');
    set(findobj('Tag','pb_zoom_out_all'),'Enable','on');
    set(findobj('Tag','pb_Autoscale'),'Enable','on');
    set(findobj('Tag','displayed'),'Enable','on');
    set(findobj('Tag', 'goto_h'), 'Enable','on');
    set(findobj('Tag', 'goto_min'), 'Enable','on'); 
    set(findobj('Tag', 'goto_sec'), 'Enable','on'); 
    set(findobj('Tag', 'edit_gotosecond'), 'Enable','on'); 
    set(findobj('Tag', 'checkbox_buttons'), 'Enable','on'); 
    set(findobj('Tag', 'goto_h'), 'String',sprintf('%02d', th_ts));
    set(findobj('Tag', 'goto_min'), 'String',sprintf('%02d', tmin_ts)); 
    set(findobj('Tag', 'goto_sec'), 'String',sprintf('%06.3f', tsec1_ts)); 
    Data.File.file = file;
    Data.File.path = path;
    set(findobj('Tag','sviewer'),'UserData',Data);
    a = Data.NS; % number of channels
    drawplot(a,only_plot);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw plots
function drawplot(numb_channel,only_plot)

    Data = get(findobj('Tag','sviewer'),'UserData');
    Data.SPlot = [];
    b = 0;
    cb_range = get(findobj('Tag','checkbox_range'),'Value');
    cb_grid = get(findobj('Tag','checkbox_grid'),'Value');
    cb_buttons = get(findobj('Tag','checkbox_buttons'),'Value');
    
    plotstart = 0.1;    % Start-Position of the plots
    f = 0.789;  % f is the vertical space requirement for the plots (% of the window)
    width = 0.84; %0.77 % width of the plots
    d = f / numb_channel;
    
% % %     if cb_range == 1
% % %         plotstart = 0.1;
% % %     else
% % %         plotstart = 0.05;
% % %         width = width + 0.1 - 0.018; % 0.86; % 0.9
% % %     end
    
    if cb_buttons == 1
        plotstart = 0.1;
    else
        width = width + 0.03;
    end

    pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
    startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
    startchannel = round(startchannel);
    k = 0;
    
    for i = startchannel : (numb_channel + startchannel - 1)
        k = k + 1;
        Data.SPlot(:,k) = subplot('Position',[plotstart 0.99-d width f/numb_channel]);
        d = d + f / numb_channel;
        b = Data.Channel{i,2};
        Data.actButton = k;
        plot_chan = plot(Data.signal(:,b),'b');
        if isfield(Data,'Detection')
            if isfield(Data.Detection,'Start')
                if isequal(Data.Detection.Start,'on')
                    set(plot_chan,'ButtonDownFcn','sviewer(''Drawline_Callback'',gcbo,[],guidata(gcbo))','UserData', Data.actButton);
                end
            end
        end
        plot_line(b);
        set_yTick(Data,b,numb_channel,cb_range);
        %set_yTick(Data,k,numb_channel,cb_range);
        if cb_grid == 1
            grid on;
        end 

        anz_xticks = length(get(gca, 'xtick'));
        if k == numb_channel
            try
                actpossl = Data.Slider.tsec;
            catch
                actpossl = 0;
            end
			x_ticklabel1 = [actpossl:Data.NoS/10:Data.NoS+actpossl];
			x_ticklabel = round(x_ticklabel1 * 100) / 100;
			
			if size(x_ticklabel,2) > size(unique(x_ticklabel),2)
                min_xtick = abs(min(x_ticklabel1));
                if min_xtick == 0
                    x_temp=x_ticklabel1(find(x_ticklabel1 > 0));
                    min_xtick = abs(min(x_temp));
                end
                factor = abs(ceil(log10(min_xtick-floor(min_xtick)))) + 1;
                x_ticklabel = round(x_ticklabel1*100*10^factor)/(100*10^factor);
			end
            set(gca, ...
                'UserData', Data.actButton, ...
                'Units', 'Normalized', ...
                'XLim',[0 Data.NoS*max(Data.HDR.SampleRate)], ...
                'XTick',[0:Data.NoS*max(Data.HDR.SampleRate)/10:Data.NoS*max(Data.HDR.SampleRate)], ...
                'XTickLabel',x_ticklabel, ...
                'YLim', [Data.ChannelConf.Display_min(b) Data.ChannelConf.Display_max(b)]);

            if isfield(Data,'Detection')
                if isfield(Data.Detection,'Start')
                    if isequal(Data.Detection.Start,'on')
                        set(gca,'ButtonDownFcn','sviewer(''Drawline_Callback'',gcbo,[],guidata(gcbo))');
                    end
                end
            end
        else
            set(gca, ...
                'UserData', Data.actButton, ...
                'XLim',[0 Data.NoS*max(Data.HDR.SampleRate)], ...
                'XTick',[0:Data.NoS*max(Data.HDR.SampleRate)/10:Data.NoS*max(Data.HDR.SampleRate)], ...
                'XTickLabel','', ...
                'YLim', [Data.ChannelConf.Display_min(b) Data.ChannelConf.Display_max(b)]); %, ...
            if isfield(Data,'Detection')
                if isfield(Data.Detection,'Start')
                    if isequal(Data.Detection.Start,'on')
                        set(gca,'ButtonDownFcn','sviewer(''Drawline_Callback'',gcbo,[],guidata(gcbo))');
                    end
                end
            end
        end
        
        if only_plot == 0
            if size(Data.HDR.PhysDim) > 1
                text_dim = Data.HDR.PhysDim{b};
            else
                text_dim = Data.HDR.PhysDim{1};
            end  
            text_dim = deblank(text_dim);
            
            x_koord = Data.NoS*max(Data.HDR.SampleRate)+Data.NoS*max(Data.HDR.SampleRate)*0.03;
            plotpos = get(gca,'Position');
            y_koord = plotpos(2);
	
            text_label = deblank(Data.Channel{i,1});
            if cb_range == 1
                text_dim = strvcat(text_label,text_dim);
                if Data.NS > 8
                    ylabel(text_label,'Rotation',0,'HorizontalAlignment', 'right','Color',[0, 0.459, 0]);
                else
                    ylabel(text_dim,'Color',[0, 0.459, 0]);
                end
            else
                if Data.NS > 8
                    ylabel(text_label,'Rotation',0,'HorizontalAlignment', 'right','Color',[0, 0.459, 0]);
                else
                    ylabel(text_label,'Color',[0, 0.459, 0]);
                end
            end
            
            x1 = 0.942; %0.933; % Start-Position of the buttons
            x2 = y_koord + (f/numb_channel)/2;
            x3 = 0.013;
            x4 = ((f-0.05)/numb_channel)/2; %0.015;
            
            if rem(k,2) == 0
                bcolor = [0.502, 0.502, 1];
            else
                bcolor = [0, 0.502, 0];
            end
            
            if cb_buttons == 1
                bcolor = [0, 0.502, 1]; % Button-Color
                fcolor = [1, 1, 0]; % Button-Font Color
                fsize = 18;
                Data.Label.Background(:,k) = uicontrol(gcf, ...
                    'Style', 'Frame', ...
                    'Units', 'Normalized', ...
                    'Position', [x1-0.004,x2-((f-0.39)/numb_channel),x3*2+0.01,((f-0.39)/numb_channel)*2], ...
                    'BackgroundColor', [1,1,1], ...
                    'UserData', Data.actButton);
                Data.Label.Zoom_in(:,k) = uicontrol(gcf, ...
                    'Style', 'PushButton', ...
                    'Units', 'normalized', ...
                    'Position', [x1,x2,x3,x4], ...
                    'String', '+', ...
                    'ForegroundColor', fcolor, ...
                    'FontWeight', 'Bold', ...
                    'FontUnit', 'pixel', ...
                    'FontName', 'Arial', ...
                    'FontSize', fsize, ...
                    'BackgroundColor', bcolor, ...
                    'Callback', 'sviewer(''LocalRescaleIN_Callback'',gcbo,[],guidata(gcbo))', ...
                    'UserData', Data.actButton);
                Data.Label.Zoom_out(:,k) = uicontrol(gcf, ...
                    'Style', 'PushButton', ...
                    'Units', 'Normalized', ...
                    'Position', [x1,x2-x4,x3,x4], ...
                    'String', '-', ...
                    'ForegroundColor', fcolor, ...
                    'FontWeight', 'Bold', ...
                    'FontUnit', 'pixel', ...
                    'FontName', 'Arial', ...
                    'FontSize', fsize, ...
                    'BackgroundColor', bcolor, ...
                    'Callback', 'sviewer(''LocalRescaleOUT_Callback'',gcbo,[],guidata(gcbo))', ...
                    'UserData', Data.actButton);
                Data.Label.Info(:,k) = uicontrol(gcf, ...
                    'Style', 'PushButton', ...
                    'Units', 'Normalized', ...
                    'Position', [x1+x3,x2-x4,x3,x4*2], ...
                    'String', '?', ...
                    'ForegroundColor', fcolor, ...
                    'FontWeight', 'Bold', ...
                    'FontUnit', 'pixel', ...
                    'FontName', 'Arial', ...
                    'FontSize', fsize, ...
                    'BackgroundColor', bcolor, ...
                    'Callback', 'sviewer(''Channel_conf_Callback'',gcbo,[],guidata(gcbo))', ...
                    'UserData', Data.actButton);
            end
            try delete(Data.Patch.firstline); end
            Data.Patch = [];
        end
    end
set(findobj('Tag','sviewer'),'UserData',Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate YTICK
function set_yTick(Data,b,numb_channel,cb_range)

akt_min = Data.ChannelConf.Display_min(b);
akt_max = Data.ChannelConf.Display_max(b);

if numb_channel <= 8
    basisline = (akt_max + akt_min)/2;
    step = (akt_max - akt_min)/8;
    y_lim = [akt_min akt_max];
    y_tick1 = [akt_min:step:akt_max];
    y_tick = round(y_tick1(2:end-1)*100)/100;
    if size(y_tick,2) > size(unique(y_tick),2)
        min_ytick = abs(min(y_tick1));
        if min_ytick == 0
            y_temp = y_tick1(find(y_tick1 > 0));
            min_ytick = abs(min(y_temp));
        end
        factor = abs(ceil(log10(step-floor(step)))) + 1;
        y_tick = round(y_tick1(2:end-1)*100*10^factor)/(100*10^factor);
    end
    y_ticklabel = y_tick;
    set(findobj('Tag','checkbox_range'),'Enable','on');
end

if numb_channel > 8
    y_lim = [akt_min akt_max];
    y_tick = y_lim;
    y_ticklabel = '';
    set(findobj('Tag','checkbox_range'),'Enable','off');
end

if size(Data.HDR.PhysDim) > 1
    text_dim = Data.HDR.PhysDim{b};
else
    text_dim = Data.HDR.PhysDim{1};
end  
text_dim = deblank(text_dim);
text_label = deblank(Data.Channel{b,1});

if cb_range == 1
    set(gca, ...
        'YLim', y_lim, ...
        'YTick', [y_tick], ...
        'YTickLabel', [y_ticklabel]);
    text_dim = strvcat(text_label,text_dim);
    if Data.NS > 8
        ylabel(text_label,'Rotation',0,'HorizontalAlignment', 'right','Color',[0, 0.459, 0]);
    else
        ylabel(text_dim,'Color',[0, 0.459, 0]);
    end
else
    set(gca, ...
        'YLim', y_lim, ...
        'YTick', [y_tick], ...
        'YTickLabel', '');
    if Data.NS > 8
        ylabel(text_label,'Rotation',0,'HorizontalAlignment', 'right','Color',[0, 0.459, 0]);
    else
        ylabel(text_label,'Color',[0, 0.459, 0]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw plot
function plot_line(channel_numb)

Data = get(findobj('Tag','sviewer'),'UserData');

if isfield(Data,'Detection')
    slider_step = get(findobj('Tag','Slider1'),'SliderStep');
    pos = Data.Slider.Pos;
    % X = startpoint of the plot in samples
    X = pos/slider_step(1)*max(Data.HDR.SampleRate)*Data.NoS/2;
    % Y = endpoint of the plot in samples
    Y = X + max(Data.HDR.SampleRate)*Data.NoS;
    if isequal(Data.Detection.Start,'on')
        if size(Data.Detection.EventMatrix) == 0
            return
        end
        found_det = find((Data.Detection.EventMatrix(:,3) == channel_numb | Data.Detection.EventMatrix(:,3) == 0 ) & ...
                    (Data.Detection.EventMatrix(:,1)<=Y) & ...
                    (Data.Detection.EventMatrix(:,1) + Data.Detection.EventMatrix(:,4) >= X));
        for i=1:size(found_det,1)
            Data.actButton=1;
            det_typ = Data.Detection.EventMatrix(found_det(i),2);
            startpoint = max(X,Data.Detection.EventMatrix(found_det(i),1));
            endpoint = min(Y,(Data.Detection.EventMatrix(found_det(i),4)));
            color=get_color(det_typ,Data);
            if ~isequal(color,[])
                id = Data.Detection.EventMatrix(found_det(i),5);
                if startpoint-X == 0
                    endpoint = abs(endpoint + Data.Detection.EventMatrix(found_det(i),1) - X);
                end
                detpatch('position',[startpoint-X -10000 endpoint 20000], ...
                        'facecolor',color, ...
                        'id', id, ...
                        'tag', ['detpatch_',int2str(id)], ...
                        'data', Data);
            end 
        end
    end
end
hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color for the patch ist defined
function color = get_color(det_typ,Data)

s = Data.Eventcodes_txt;
dettyp_dec = dec2hex(det_typ);
pos = hex2dec(dettyp_dec(end));

if pos == 0
    pos = 17;
end
select = Data.ShowDetection.MainClass;
select_string = Data.ShowDetection.SubClass;

if strmatch('all Events',select_string,'exact')
    if (select + 1) > length(s.GroupValue)
        if (det_typ >= s.GroupValue(select) & det_typ < 32768)
            color = Data.detcolor.detcolor(pos,:);
        else
            color = [];
        end
    else
        if s.GroupValue(select) == 1056
            if (det_typ >= s.GroupValue(select) & det_typ < s.GroupValue(select + 1) | det_typ == 33824)
                color = Data.detcolor.detcolor(pos,:);
            else
                color = [];
            end
        else
            if (det_typ >= s.GroupValue(select) & det_typ < s.GroupValue(select + 1))
                color = Data.detcolor.detcolor(pos,:);
            else
                color = [];
            end
        end
    end
else
    pos_desc = find(s.CodeIndex == det_typ);
    if strmatch(s.CodeDesc(pos_desc),select_string,'exact')
        color = Data.detcolor.detcolor(pos,:);
    else
        color = [];
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare four intervalls
function [startpos,endpos,find] = comp_intervall(K,L,X,Y)
if K>=X & K<=Y & L>=X & L<=Y
    startpos = K;
    endpos = L;
    find=true;
    return;
end
if K<X & L>X & L<=Y
    startpos = X;
    endpos = L;
    find=true;
    return;
end
if K>=X & K<Y & L>Y
    startpos = K;
    endpos = Y;
    find=true;
    return;
end
if K<X & L>Y
    startpos = X;
    endpos = Y;
    find=true;
    return;
end
if K>Y & L>Y
    startpos = -1;
    endpos = -1;
    find = false;
else
    startpos = -1;
    endpos = -1;
    find=true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete old Buttons, Text....
function deleteObj(Data,only_plot)
if ~isempty(Data);
    if only_plot == 0
        try delete(Data.Label.Zoom_out); end
        try delete(Data.Label.Zoom_in); end
        try delete(Data.Label.Background); end
        try delete(Data.Label.Info); end
        try delete(Data.SPlot); end
        try Data.Label.Text=[]; end
        try Data.Label.Zoom_out=[]; end
        try Data.Label.Zoom_in=[]; end
        try Data.Label.Info=[]; end
        try Data.Label.Background=[]; end
        try Data.SPlot=[]; end
        try Data.Patch = []; end
        if isfield(Data.Label,'Range')
            delete(Data.Label.Range);
            Data.Label.Range=[];
        end
    else
        try delete(Data.SPlot); end
        try Data.SPlot=[]; end
    end
    set(findobj('Tag','sviewer'),'UserData',Data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% redraw plots after changing
function redraw_plots(h,numb_plot,Data,aktpos)

pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
startchannel = round(startchannel);

for i = 1:numb_plot
    subplot(h(i));
    set(gca,'NextPlot','replacechildren');
    b = Data.Channel{startchannel,2};
    startchannel = startchannel + 1;
    plot_chan = plot(Data.signal(:,b),'b');
    
    cb_range = get(findobj('Tag','checkbox_range'),'Value');
    set_yTick(Data,b,numb_plot,cb_range)
    %set_yTick(Data,i,numb_plot,cb_range);
    
    set(findobj('Tag','sviewer'),'UserData',Data);

    if i == numb_plot
		x_ticklabel1 = [aktpos:Data.NoS/10:Data.NoS+aktpos];
		x_ticklabel = round(x_ticklabel1 * 100) / 100;
		if size(x_ticklabel,2) > size(unique(x_ticklabel),2)
            min_xtick = abs(min(x_ticklabel1));
            if min_xtick == 0
                x_temp=x_ticklabel1(find(x_ticklabel1 > 0));
                min_xtick = abs(min(x_temp));
            end
            factor = abs(ceil(log10(min_xtick-floor(min_xtick)))) + 1;
            x_ticklabel = round(x_ticklabel1*100*10^factor)/(100*10^factor);
		end
        set(gca, ...
            'XTick',[0:Data.NoS*max(Data.HDR.SampleRate)/10:Data.NoS*max(Data.HDR.SampleRate)], ...
            'XTickLabel',x_ticklabel);
    end
    if isfield(Data,'Detection')
        if isfield(Data.Detection,'Start')
            if isequal(Data.Detection.Start,'on')
                set(plot_chan,'ButtonDownFcn','sviewer(''Drawline_Callback'',gcbo,[],guidata(gcbo))','UserData', Data.actButton);
            end
        end
    end
    plot_line(b);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function Slider1_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag','sviewer'),'UserData');
length = Data.Total_length_samples;
last_pos = Data.Slider.Pos;
pos1 = get(hObject, 'Value');

numb_samples = Data.Total_length_samples;
max_length = numb_samples - Data.ShowSamples;
step = Data.ShowSamples / max_length / 2; 

if rem(pos1,step) > 0
    new_pos = round(pos1/step);
    pos1 = min(new_pos*step,1);
    set(findobj('Tag', 'Slider1'), 'Value',pos1); 
end

pos = pos1;
if pos <= (length-Data.ShowSamples) 
    if pos ~= last_pos
        tsec = pos * (length-Data.ShowSamples) / max(Data.HDR.SampleRate);
        tmin = floor(tsec/60);
        tsec1 = rem(tsec,60);
        th = floor(tmin/60);
        tmin = rem(tmin,60);
        set(findobj('Tag', 'goto_h'), 'String',sprintf('%02d', th));
        set(findobj('Tag', 'goto_min'), 'String',sprintf('%02d', tmin)); 
        set(findobj('Tag', 'goto_sec'), 'String',sprintf('%06.3f', tsec1));
        read_dur = tsec + Data.NoS;
        if read_dur > Data.Total_length_sec
            tsec = Data.Total_length_sec - Data.NoS;
            tsec = round(tsec*Data.HDR.SampleRate)/Data.HDR.SampleRate;
            [Data.signal,Data.HDR] = sread(Data.HDR,Data.NoS,tsec);
        else
            tsec = round(tsec*Data.HDR.SampleRate)/Data.HDR.SampleRate;
            [Data.signal,Data.HDR] = sread(Data.HDR,Data.NoS,tsec);
        end
        % h = position of the plots
        Data.Slider.Pos = pos;
        set(findobj('Tag','sviewer'),'UserData',Data);
        h = Data.SPlot;
        numb_plot = max(size(h));
        redraw_plots(h,numb_plot,Data,tsec);
        Data.Slider.tsec = tsec;
        set(findobj('Tag','sviewer'),'UserData',Data);
    end
end

% if the startpoint of the detection is marked and in the meantime the
% slider is used the startpoint is drawn again
try
	if isfield(Data.Patch,'firstline')
        if ~isequal(Data.Patch.firstline,[])
            slider_step = get(findobj('Tag','Slider1'),'SliderStep');
            pos = Data.Slider.Pos;
            X = pos/slider_step(1)*max(Data.HDR.SampleRate)*Data.NoS / 2;
            pos = Data.Patch.line;
            if X < pos
                h = Data.SPlot;
                pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
                startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
                startchannel = round(startchannel);
                subplot(h(Data.Patch.handleline-startchannel+1));
                Data.Patch.firstline = line('XData', [pos-X pos-X], ...
                    'YData', [-99999 99999], ...
                    'EraseMode', 'xor', ...
                    'Color', [1 0.2 0.2]);
                set(findobj('Tag','sviewer'),'UserData',Data);
            end
        end
	end
end
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function Slider_Channel_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');

if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end
file = Data.File.file;
path = Data.File.path;
drawnew (file,path,Data.ShowSamples,1);
Data = get(findobj('Tag', 'sviewer'), 'UserData');
% if the startpoint of the detection is marked and in the meantime the
% slider is used the startpoint is drawn again
try
	if isfield(Data.Patch,'firstline')
        if ~isequal(Data.Patch.firstline,[])
            slider_step = get(findobj('Tag','Slider1'),'SliderStep');
            
            pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
            startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
            startchannel = round(startchannel);
            
            pos = Data.Slider.Pos;
            X = pos/slider_step(1)*max(Data.HDR.SampleRate)*Data.NoS / 2;
            pos = Data.Patch.line;
            if (X < pos) & (Data.Patch.selChannel >= startchannel) & (Data.Patch.selChannel <= startchannel + Data.NS -1)
                h = Data.SPlot;
                subplot(h(Data.Patch.handleline-startchannel+1));
                Data.Patch.firstline = line('XData', [pos-X pos-X], ...
                    'YData', [-99999 99999], ...
                    'EraseMode', 'xor', ...
                    'Color', [1 0.2 0.2]);
                set(findobj('Tag','sviewer'),'UserData',Data);
            end
        end
	end
end

figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Slider_ShowChannels_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');

pos1 = get(hObject, 'Value');
numb = pos1 * size(Data.Channel,1);

if numb > 0 & numb ~= str2num(get(findobj('Tag','numb_act_channels'),'String'))
    set(findobj('Tag','numb_act_channels'),'String',numb);
    file = Data.File.file;
    path = Data.File.path;
    try delete(Data.Patch.firstline); end
    Data.Patch = [];
    set(findobj('Tag','sviewer'),'UserData',Data);
    drawnew (file,path,Data.ShowSamples,0);
end

figure(gcf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pb_Autoscale_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');

if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end

for k = 1:length(Data.ChannelConf.Display_min)
    sample_min = min(Data.signal(:,k));
    sample_max = max(Data.signal(:,k));
    if sample_min == sample_max
        sample_min = sample_min -10;
        sample_max = sample_max +10;
    end
    Data.ChannelConf.Display_min(k,:) = sample_min;
    Data.ChannelConf.Display_max(k,:) = sample_max;
end

set(findobj('Tag','sviewer'),'UserData',Data);
h = Data.SPlot;
numb_plot = max(size(h));
redraw_plots(h,numb_plot,Data,Data.Slider.tsec);
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Edit_gotosecond_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');

maxtime = Data.Total_length_sec;
gotime = get(hObject,'String');

if isempty(gotime)
    errordlg('No value for Go To Second!', 'Error');
    return;
end
gotime = str2double(gotime);
if gotime > maxtime
    errordlg(['The value for Go To Second is too large ! [max. ' int2str(maxtime) ' sec.]'], 'Error');
    set(hObject,'String','0');
    return;
end
if isnan(gotime)
    errordlg('The value for Go To Second is Not a Number!', 'Error');
    set(hObject,'String','0');
    return;
end
if gotime >= (maxtime-Data.NoS)
    gotime = gotime - Data.NoS;
end
Data.Display.gotime = gotime;
set(findobj('Tag','sviewer'),'UserData',Data);

file = Data.File.file;
path = Data.File.path;
length = Data.Total_length_sec;

pos = gotime/(length-Data.NoS);
set(findobj('Tag', 'Slider1'), 'Value',pos); 
last_pos = Data.Slider.Pos;
if pos <= (length-Data.NoS) 
    if pos ~= last_pos
        tsec = round(pos*(length-Data.NoS));
        tmin = floor(tsec/60);
        tsec1 = rem(tsec,60);
        th = floor(tmin/60);
        tmin = rem(tmin,60); 
        set(findobj('Tag', 'goto_h'), 'String',sprintf('%02d', th));
        set(findobj('Tag', 'goto_min'), 'String',sprintf('%02d', tmin)); 
        set(findobj('Tag', 'goto_sec'), 'String',sprintf('%06.3f', tsec1));
        [Data.signal,Data.HDR] = sread(Data.HDR,Data.NoS,gotime);
        h = Data.SPlot;
        numb_plot = max(size(h));
        Data.Slider.Pos = pos;
        set(findobj('Tag','sviewer'),'UserData',Data);
        redraw_plots(h,numb_plot,Data,tsec);
        Data.Slider.Pos = pos;
        Data.Slider.tsec = tsec;
        set(findobj('Tag','sviewer'),'UserData',Data);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change the zoom-factor
function Zoom_factor_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
scale_factor = str2num(get(hObject,'String')); 
if length(scale_factor) > 1 | length(scale_factor) == 0
    errordlg('The value for Zoom-Factor is Not a Number!', 'Error');
    set(hObject,'String','2');
    return;
end
if scale_factor <= 1
    errordlg('Zoom-Factor must be larger than 1 (e.g. 1.1) !', 'Error');
    set(hObject,'String','2');
    return;
end
Data.ChannelConf.Scale = str2num(get(hObject,'String')); 
set(findobj('Tag','sviewer'),'UserData',Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Checkbox_buttons_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
file = Data.File.file;
path = Data.File.path;
drawnew (file,path,Data.ShowSamples,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to change the diplay-time or jump to a indicated second of the signal
function DisplayMenu_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end

waitfor(sviewer_display);
Data = get(findobj('Tag', 'sviewer'), 'UserData');
try delete(Data.Patch.firstline); end
Data.Patch = [];
set(findobj('Tag','sviewer'),'UserData',Data);
displaytrue = Data.Display;

if isempty(displaytrue)
    return;
end

file = Data.File.file;
path = Data.File.path;
length = Data.Total_length_sec;
displaytime = Data.Display.newdisplaytime;
gotime = Data.Display.gotime;

channel=Data.Channel;
if isempty(channel)
    channel = {};
end

if ~isempty(displaytime)
    set(findobj('Tag', 'Slider1'), 'Value',0);
    Data.ShowSamples = displaytime * max(Data.HDR.SampleRate);
    set(findobj('Tag','sviewer'),'UserData',Data);
    drawnew (file,path,Data.ShowSamples,0);
end

if ~isempty(gotime)
    pos=gotime/(length-Data.NoS);
    set(findobj('Tag', 'Slider1'), 'Value',pos); 
    last_pos = Data.Slider.Pos;
    if pos <= (length-Data.NoS) 
        if pos ~= last_pos
            tsec = round(pos*(length-Data.NoS));
            tmin = floor(tsec/60);
            tsec1 = rem(tsec,60);
            th = floor(tmin/60);
            tmin = rem(tmin,60); 
            set(findobj('Tag', 'goto_h'), 'String',sprintf('%02d', th));
            set(findobj('Tag', 'goto_min'), 'String',sprintf('%02d', tmin)); 
            set(findobj('Tag', 'goto_sec'), 'String',sprintf('%06.3f', tsec1));
            [Data.signal,Data.HDR] = sread(Data.HDR,Data.NoS,gotime);
            h = Data.SPlot;
            numb_plot = max(size(h));
            Data.Slider.Pos = pos;
            set(findobj('Tag','sviewer'),'UserData',Data);
            redraw_plots(h,numb_plot,Data,tsec);
            Data.Slider.Pos = pos;
            Data.Slider.tsec = tsec;
            set(findobj('Tag','sviewer'),'UserData',Data);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose the channels you want to display
function ChannelsMenu_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');

if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end

waitfor(sviewer_channel);
Data = get(findobj('Tag', 'sviewer'), 'UserData');
try delete(Data.Patch.firstline); end
Data.Patch = [];
set(findobj('Tag','sviewer'),'UserData',Data);
newchannel = Data.Channel;

% % % if isequal(Data.changeChannel,'no')
% % %     return;
% % % end

set(findobj('Tag','sviewer'),'UserData',Data);
file = Data.File.file;
path = Data.File.path;

% % % size_newchannel = size(newchannel);
% % % size_allchannel = size(Data.allChannel);
% % % 
% % % for i = 1:size_newchannel(1)
% % %     sel_channel = deblank(Data.Channel(i,:));
% % %     b = 0;
% % %     search = true;
% % %     while search
% % %         b = b + 1;
% % %         data_channel = deblank(Data.allChannel(b,:));
% % %         if isequal(sel_channel,data_channel)
% % %             search = false;
% % %             text = data_channel{1,1};
% % %         end
% % %     end
% % %     channel{i,1}=text;
% % %     channel{i,2}=b;
% % % end
% % % Data.Channel = {};
% % % Data.Channel = channel;
% % % 
% % % set(findobj('Tag','numb_act_channels'),'String',min(Data.ShowChannelmax,size_newchannel(1)));

set(findobj('Tag','sviewer'),'UserData',Data);
drawnew (file,path,Data.ShowSamples,0);
all_channels{1}='all Channels';
show_channels={};
show_channels=[all_channels;Data.Channel(:,1)];
set(findobj('Tag','listbox_detchannels'),'Value',1);
set(findobj('Tag','listbox_detchannels'),'String',show_channels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Edit_Displayed_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');

if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end

try delete(Data.Patch.firstline); end
Data.Patch = [];
set(findobj('Tag','sviewer'),'UserData',Data);
Data.Display.newdisplaytime = '';
maxtime = Data.Total_length_sec;
file = Data.File.file;
path = Data.File.path;

displaytime = str2double(get(findobj('Tag', 'displayed'), 'String'));
if displaytime > maxtime
    errordlg(['The value for Displayed-Time is too large ! [max. ' num2str(maxtime) ' sec.]'], 'Error');
    set(findobj('Tag', 'displayed'), 'String',num2str(Data.NoS));
    return;
end
if isnan(displaytime)
    errordlg('The value for Displayed-Time is Not a Number!', 'Error');
    set(findobj('Tag', 'displayed'), 'String',num2str(Data.NoS));
    return;
end
if length(displaytime) > 1 | length(displaytime) == 0
    errordlg('The value for Display-Time is Not a Number!', 'Error');
    set(findobj('Tag', 'displayed'), 'String',num2str(Data.NoS));
    return;
end
if displaytime == 0
    errordlg('The value must be larger than 0!', 'Error');
    set(findobj('Tag', 'displayed'), 'String',num2str(Data.NoS));
    return;
end
Data.Display.newdisplaytime = displaytime;
channel = Data.Channel;
if isempty(channel)
    channel = {};
end

if ~isempty(displaytime)
    Data.ShowSamples = displaytime * max(Data.HDR.SampleRate);
    set(findobj('Tag','sviewer'),'UserData',Data);
    drawnew (file,path,Data.ShowSamples,1);
end
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Edit_Goto_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');

if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end

try delete(Data.Patch.firstline); end
Data.Patch = [];
set(findobj('Tag','sviewer'),'UserData',Data);
Data.Display.gotime = '';
file = Data.File.file;
path = Data.File.path;
length = Data.Total_length_sec;

goto_h = str2num(get(findobj('Tag', 'goto_h'), 'String'));
goto_min = str2num(get(findobj('Tag', 'goto_min'), 'String'));
goto_sec = str2num(get(findobj('Tag', 'goto_sec'), 'String'));

pos_s1 = get(findobj('Tag', 'Slider1'), 'Value');
Data.Slider.Pos = pos_s1;
length_s1 = Data.Total_length_samples;
tsec_s1 = pos_s1 * (length_s1 - Data.ShowSamples) / max(Data.HDR.SampleRate);
tmin_ts = floor(tsec_s1 / 60);
tsec1_ts = rem(tsec_s1,60);
th_ts = floor(tmin_ts / 60);
tmin_ts = rem(tmin_ts,60); 
tsec_s1 = round(tsec_s1 * Data.HDR.SampleRate) / Data.HDR.SampleRate;

if size(goto_h,2) > 1 | size(goto_min,2) > 1 | size(goto_sec,2) > 1
    errordlg('The value for Go To Second is wrong!', 'Error');
    return;
end

if isempty(goto_h) | isempty(goto_min) | isempty(goto_sec)
    errordlg('The value for Go To Second is wrong!', 'Error');
    return;
end

gotime = (goto_h * 3600) +  (goto_min * 60) + goto_sec;

if isempty(gotime)
    errordlg('No value for Go To Second!', 'Error');
    return;
end
if gotime > length
    tsec = length;
    tmin = floor(tsec/60);
    tsec1 = rem(tsec,60);
    th = floor(tmin/60);
    tmin = rem(tmin,60);
    max_value = sprintf('%02d:%02d:%06.3f', th,tmin,tsec1);
    errordlg(['The value for Go To Second is too large ! [max. ' max_value ' sec.]'], 'Error');
    set(findobj('Tag', 'goto_h'), 'String',sprintf('%02d', th_ts));
    set(findobj('Tag', 'goto_min'), 'String',sprintf('%02d', tmin_ts));
    set(findobj('Tag', 'goto_sec'), 'String',sprintf('%06.3f', tsec1_ts));
    return;
end
if isnan(gotime)
    errordlg('The value for Go To Second is Not a Number!', 'Error');
    return;
end
if gotime >= (length-Data.NoS)
    gotime = gotime - Data.NoS;
end
Data.Display.gotime = gotime;

channel = Data.Channel;
if isempty(channel)
    channel = {};
end

if ~isempty(gotime)
    pos = gotime/(length-Data.NoS);
    set(findobj('Tag', 'Slider1'), 'Value',pos); 
    last_pos = Data.Slider.Pos;
    if pos <= (length-Data.NoS) 
        if pos ~= last_pos
            tsec = round(pos*(length-Data.NoS));
            tmin = floor(tsec/60);
            tsec1 = rem(tsec,60);
            th = floor(tmin/60);
            tmin = rem(tmin,60); 
            set(findobj('Tag', 'goto_h'), 'String',sprintf('%02d', th));
            set(findobj('Tag', 'goto_min'), 'String',sprintf('%02d', tmin)); 
            set(findobj('Tag', 'goto_sec'), 'String',sprintf('%06.3f', tsec1));
            [Data.signal,Data.HDR] = sread(Data.HDR,Data.NoS,gotime);
            h = Data.SPlot;
            numb_plot = max(size(h));
            Data.Slider.Pos = pos;
            set(findobj('Tag','sviewer'),'UserData',Data);
            redraw_plots(h,numb_plot,Data,tsec);
            Data.Slider.Pos = pos;
            Data.Slider.tsec = tsec;
            set(findobj('Tag','sviewer'),'UserData',Data);
        end
    end
end

set(findobj('Tag','sviewer'),'UserData',Data);
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configurate the active channel
function Channel_conf_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end
try delete(Data.Patch.firstline); end
Data.Patch = [];
set(findobj('Tag','sviewer'),'UserData',Data);
Data.Channelconf.whichbutton = get(gcbo, 'UserData');
set(findobj('Tag', 'sviewer'), 'UserData',Data);
pause(0.0000000001);
waitfor(sviewer_channel_conf);
Data = get(findobj('Tag', 'sviewer'), 'UserData');
Data.Channelconf.whichbutton = [];
h = Data.SPlot;
numb_plot = max(size(h));
redraw_plots(h,numb_plot,Data,Data.Slider.tsec);
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zoom in
function LocalRescaleIN_Callback(hObject, eventdata, handles)
LocalRescale('in',0);   % 0... zoom in -> act. Channel
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zoom out
function LocalRescaleOUT_Callback(hObject, eventdata, handles)
LocalRescale('out',0);  % 0... zoom in -> act. Channel
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zoom out all channels
function pb_zoom_out_all_Callback(hObject, eventdata, handles)
LocalRescale('out',1);  % 1... zoom in -> all Channels
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zoom in all channels
function pb_zoom_in_all_Callback(hObject, eventdata, handles)
LocalRescale('in',1);  % 1... zoom in -> all Channels
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the new values (Zoom-in / Zoom-out)
function LocalRescale(Mode,settings)

% Mode: in / out
% settings: 0 ... act. Channel
%           1 ... all Channels

Data = get(findobj('Tag','sviewer'),'UserData');

pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
startchannel = round(startchannel);

if settings == 0
    whichbutton = get(gcbo, 'UserData');
    actchannel = whichbutton + startchannel - 1;
    numb_channels = 1;
else
    actchannel = startchannel;
    numb_channels = Data.NS;
end

length = Data.Total_length_sec;
last_pos = Data.Slider.Pos;
pos1 = get(findobj('Tag','Slider1'), 'Value');
max_length = length - Data.NoS;
step = Data.NoS/max_length;
if rem(pos1,step) > 0
    new_pos = round(pos1/step);
    pos1 = new_pos*step;
end
pos = pos1;
aktpos = round(pos*(length-Data.NoS));

for i = 1:numb_channels
	b = Data.Channel{actchannel,2};
	disp_min = Data.ChannelConf.Display_min(b);
	disp_max = Data.ChannelConf.Display_max(b);
	scale_factor = Data.ChannelConf.Scale;
	basisline = (disp_max + disp_min)/2;
	switch Mode,
        case 'in',
            Data.ChannelConf.Display_min(b) = (disp_min - basisline) / scale_factor + basisline;
            Data.ChannelConf.Display_max(b) = (disp_max - basisline) / scale_factor + basisline;
        case 'out',
            Data.ChannelConf.Display_min(b) = (disp_min - basisline) * scale_factor + basisline;
            Data.ChannelConf.Display_max(b) = (disp_max - basisline) * scale_factor + basisline;
	end
	set(findobj('Tag','sviewer'),'UserData',Data);
	h = Data.SPlot;
	numb_plots = size(h);
    try
        subplot(h(whichbutton));
    catch
        subplot(h(i));
    end
    set(gca,'YLim', [disp_min disp_max]);
	cb_range = get(findobj('Tag','checkbox_range'),'Value');
	set_yTick(Data,b,numb_plots(2),cb_range);
    actchannel = actchannel + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diplay the range
function Checkbox_Range_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');

active = get(hObject,'Value');

pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
startchannel = round(startchannel)-1;

h = Data.SPlot;
numb_plots = size(h,2);
for i = 1:numb_plots
    subplot(h(i));
    set_yTick(Data,i+startchannel,numb_plots,active)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the grid
function Checkbox_Grid_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
cb_grid = get(findobj('Tag','checkbox_grid'),'Value');
h = Data.SPlot;
numb_plots = size(h);
if cb_grid == 1
    for i = 1:numb_plots(2)
        subplot(h(i));
        grid on;
    end
else
    for i = 1:numb_plots(2)
        subplot(h(i));
        grid off;
    end
end
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of displayed channels
function ShowChannels_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
try delete(Data.Patch.firstline); end
Data.Patch = [];
set(findobj('Tag','sviewer'),'UserData',Data);
numb_showchannels = str2num(get(hObject,'String'));
len_channels = size(Data.Channel,1);
if numb_showchannels > size(Data.Channel,1) | numb_showchannels <= 0
    errordlg(['The value must be between 0 and ' num2str(len_channels + 1) ' !'], 'Error');
    set(hObject,'String',num2str(length(Data.SPlot)));
    return;
else
    Data.ShowChannelmax = numb_showchannels;
    file = Data.File.file;
    path = Data.File.path;
    drawnew (file,path,Data.ShowSamples,0);
    figure(gcf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close sviewer
function Exit_Callback(hObject, eventdata, handles)

selection = questdlg(['Close ' get(handles.sviewer,'Name') '?'],...
                     ['Close ' get(handles.sviewer,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end
delete(handles.sviewer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open the File Info Window
function FileInfoMenu_Callback(hObject, eventdata, handles)

Data=get(findobj('Tag','sviewer'),'UserData');
if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end
sviewer_fileinfo(Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Startdetection_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag','sviewer'),'UserData');
if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end                  
if isequal(get(findobj('Tag','Startdetection'),'Label'),'Start                               F5')
    Event_string = sprintf('%s |',Data.Eventcodes_txt.GroupDesc{:});
    Event_string = Event_string(1:end-2);
    Data.Detection.Event_string = Event_string;
    
    if isfield(Data.Detection,'EventMatrix')
        typ = Data.Detection.EventMatrix(1,2);
    else
        typ = hex2dec('0101');    
    end
    s = Data.Eventcodes_txt;
    select = find(s.GroupValue == bitand(s.GroupMask,typ));
    pos_eventdetail = find(bitand(s.CodeIndex,s.GroupMask(select))==s.GroupValue(select));        % decoding of eventtable 

    Eventdetail_string = sprintf('%s |',Data.Eventcodes_txt.CodeDesc{pos_eventdetail});
    Eventdetail_string = Eventdetail_string(1:end-2);
    Data.Detection.Eventdetail_string = Eventdetail_string;
    set(findobj('Tag','sviewer'),'UserData',Data);
    try
        setdefault_detection_start(Data,select);
    catch
        setdefault_detection_start(Data,1);
    end
else
    setdefault_detection_stop(Data,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting default values for starting detection
function setdefault_detection_start(Data,select)

set(findobj('Tag','listbox_Event'),'Value',select);
set(findobj('Tag','Startdetection'),'Label','Stop/Save Detection      F5')
set(findobj('Tag','Event'),'Enable','on');
set(findobj('Tag','Event_detail'),'Enable','on');
set(findobj('Tag','listbox_detchannels'),'Enable','on');
set(findobj('Tag','listbox_Event'),'Enable','on');
set(findobj('Tag','listbox_Event_detail'),'Enable','on');
set(findobj('Tag', 'Detection_radiobutton_detection'), 'Enable', 'on');
set(findobj('Tag', 'Detection_radiobutton_display'), 'Enable', 'on');
set(findobj('Tag','Button_Update'),'Enable','on');
    
all_channels{1}='all Channels';
show_channels={};
show_channels=[all_channels;Data.Channel(:,1)];
set(findobj('Tag','listbox_detchannels'),'Value',1);
set(findobj('Tag','listbox_detchannels'),'String',show_channels);

Data.Detection.Start='on';
if ~isfield(Data.Detection,'EventMatrix')
    Data.Detection.EventMatrix = [];
end
Data.ShowDetection.MainClass = 1;
Data.ShowDetection.SubClass = 'all Events';
set(findobj('Tag', 'Detection_radiobutton_detection'), 'Value', 0);
set(findobj('Tag', 'Detection_radiobutton_display'), 'Value', 1);
Data.Slider.tsec=0;

set(findobj('Tag', 'Slider1'), 'Value',0);
set(findobj('Tag','sviewer'),'UserData',Data);

file = Data.File.file;
path = Data.File.path;
drawnew (file,path,Data.ShowSamples,0);
Radiobutton_Display_Callback;
Radiobutton_Detection_Callback;
set(findobj('Tag', 'Detection_radiobutton_detection'), 'Value', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting default values for stoping detection
function setdefault_detection_stop(Data,pwarn)

eventmatrix=Data.Detection.EventMatrix;
if ~isempty(eventmatrix)
    if pwarn == 1
        button = questdlg('Save Detection?','Warning');
    else
        button = 'Yes';
    end
    if isequal(button,'Cancel')
        return;
    end
    if isequal(button,'No')
        set(findobj('Tag','Startdetection'),'Label','Start                               F5');
        Data.Detection.EventMatrix = [];
        set(findobj('Tag','sviewer'),'UserData',Data);
        h = Data.SPlot;
        numb_plot = max(size(h));
        redraw_plots(h,numb_plot,Data,Data.Slider.tsec);
        redraw_display;
        return;
    end
    if isequal(button,'Yes')
        Data.Detection.Start='off';
        filename = Data.File.file(1:end-4);
        [Y, M, D, h, m, s] = datevec(now);
        [file,path]=uiputfile(sprintf('%s-%04i%02i%02i.evt',filename,Y,M,D),'Save Detection');
        if isequal(file,0) | isequal(path,0)
            set(findobj('Tag','Startdetection'),'Label','Start                               F5');
            Data.Detection.EventMatrix = [];
            set(findobj('Tag','sviewer'),'UserData',Data);
            h = Data.SPlot;
            numb_plot = max(size(h));
            redraw_plots(h,numb_plot,Data,Data.Slider.tsec);
            redraw_display;
            return;
        else
            new_detection = [path file];

            % save EVENTS in GDF format
            H.TYPE = 'EVENT';
            [p,f,e]=fileparts(new_detection);
            H.FileName = fullfile(p,[f,'.evt']);
            H.T0 = Data.HDR.T0;
            H.Patient = Data.HDR.Patient;
            H.PID = ['Scoring of file: ', Data.HDR.FileName];
            H.RID = 'Generated with SVIEWER from http://biosig.sf.net';
            H.EVENT.SampleRate = Data.HDR.SampleRate; 
            H.EVENT.POS = eventmatrix(:,1);
            H.EVENT.TYP = eventmatrix(:,2);
            H.EVENT.CHN = eventmatrix(:,3);
            H.EVENT.DUR = eventmatrix(:,4);
            H.NS = 0; 
            H.PhysDimCode = [];
            H=sopen(H,'w'); H=sclose(H);

            Data.Detection.EventMatrix = [];
            set(findobj('Tag','sviewer'),'UserData',Data);
            h = Data.SPlot;
            numb_plot = max(size(h));
            redraw_plots(h,numb_plot,Data,Data.Slider.tsec);
        end
    end
end
redraw_display;
set(findobj('Tag','Startdetection'),'Label','Start                               F5');
set(findobj('Tag','sviewer'),'UserData',Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update the display
function redraw_display

set(findobj('Tag','Startdetection'),'Label','Start                               F5');
set(findobj('Tag','listbox_detchannels'),'Value',1);
set(findobj('Tag','listbox_detchannels'),'String',[]);
set(findobj('Tag','listbox_detchannels'),'Enable','off');
set(findobj('Tag','Event'),'Enable','off');
set(findobj('Tag','Event_detail'),'Enable','off');
set(findobj('Tag','listbox_Event'),'Enable','off');
set(findobj('Tag','listbox_Event_detail'),'Enable','off');
set(findobj('Tag', 'Detection_radiobutton_detection'), 'Enable', 'off');
set(findobj('Tag', 'Detection_radiobutton_display'), 'Enable', 'off');
set(findobj('Tag','Button_Update'),'Enable','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiple alternative of different Main-Class events
function Listbox_MainClass_Callback(hObject, eventdata, handles)

select = get(findobj('Tag','listbox_Event'), 'Value');
Data = get(findobj('Tag','sviewer'),'UserData');
try delete(Data.Patch.firstline); end
Data.Patch = [];
set(findobj('Tag','sviewer'),'UserData',Data);
s = Data.Eventcodes_txt;
pos_eventdetail = find(bitand(s.CodeIndex,s.GroupMask(select))==s.GroupValue(select));        % decoding of eventtable 
Eventdetail_string = sprintf('%s |',Data.Eventcodes_txt.CodeDesc{pos_eventdetail});
Eventdetail_string = Eventdetail_string(1:end-2);
active_radiobutton = get(findobj('Tag', 'Detection_radiobutton_detection'), 'Value');
if active_radiobutton
    set(findobj('Tag','listbox_Event_detail'),'String',Eventdetail_string,'Value',1);
else
    set(findobj('Tag','listbox_Event_detail'),'String',['all Events |' Eventdetail_string],'Value',1);
    select_string = get(findobj('Tag','listbox_Event_detail'),'String');
    select_value = get(findobj('Tag','listbox_Event_detail'),'Value');
    Data.ShowDetection.SubClass = select_string(select_value,:);
end

Data.ShowDetection.MainClass = select;
set(findobj('Tag','sviewer'),'UserData',Data);

h = Data.SPlot;
numb_plot = max(size(h));
redraw_plots(h,numb_plot,Data,Data.Slider.tsec);
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiple alternative of different Sub-Class events
function Listbox_SubClass_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag','sviewer'),'UserData');
try delete(Data.Patch.firstline); end
Data.Patch = [];
display_active = get(findobj('Tag', 'Detection_radiobutton_display'), 'Value');
if display_active
    Data = get(findobj('Tag','sviewer'),'UserData');
    select_string = get(findobj('Tag','listbox_Event_detail'),'String');
    select_value = get(findobj('Tag','listbox_Event_detail'),'Value');
    Data.ShowDetection.SubClass = select_string(select_value,:);
    set(findobj('Tag','sviewer'),'UserData',Data);
    h = Data.SPlot;
    numb_plot = max(size(h));
    redraw_plots(h,numb_plot,Data,Data.Slider.tsec);
end
set(findobj('Tag','sviewer'),'UserData',Data);
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The last executed detection can be taken over to further selected channels
function Listbox_Update_Channels_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
if ~isfield(Data.Detection,'Update')
    return
end
display_active = get(findobj('Tag', 'Detection_radiobutton_display'), 'Value');
if display_active
    return
end
try delete(Data.Patch.firstline); end
Data.Patch = [];
set(findobj('Tag','sviewer'),'UserData',Data);
values = get(findobj('Tag','listbox_detchannels'),'Value');
size_selected = length(values);
% values == 1 -> all Channels
if find(values == 1)
    i = size(Data.Detection.EventMatrix,1);
    if i == 0
        id = 1;
    else
        id = Data.Detection.EventMatrix(end,5)+1;
    end
    startline = Data.Detection.Update(1);
    det_typ = Data.Detection.Update(2);
    dur = Data.Detection.Update(4);
    exist_data = find(Data.Detection.EventMatrix(:,1) == startline & Data.Detection.EventMatrix(:,2) == det_typ & ...
                          Data.Detection.EventMatrix(:,3) >= 0 & Data.Detection.EventMatrix(:,4) == dur);
    if exist_data
        Data.Detection.EventMatrix(exist_data,:) = [];
        Data.Detection.EventMatrix(i,:) = [startline,det_typ,0,dur,id];
    else
        Data.Detection.EventMatrix(i+1,:) = [startline,det_typ,0,dur,id];
    end
else
    for j = 1:size_selected
        i = size(Data.Detection.EventMatrix,1);
        if i == 0
            id = 1;
        else
            id = Data.Detection.EventMatrix(end,5)+1;
        end
        startline = Data.Detection.Update(1);
        det_typ = Data.Detection.Update(2);
        sel_channel = Data.Channel{values(j)-1,2};
        dur = Data.Detection.Update(4);
        exist_data = find(Data.Detection.EventMatrix(:,1) == startline & Data.Detection.EventMatrix(:,2) == det_typ & ...
                          (Data.Detection.EventMatrix(:,3) == sel_channel | Data.Detection.EventMatrix(:,3) == 0) & ...
                          Data.Detection.EventMatrix(:,4) == dur);
        if length(exist_data) == 0
            Data.Detection.EventMatrix(i+1,:) = [startline,det_typ,sel_channel,dur,id];
        end
    end
end
set(findobj('Tag','sviewer'),'UserData',Data);

h = Data.SPlot;
numb_plot = max(size(h));
redraw_plots(h,numb_plot,Data,Data.Slider.tsec);
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% load events
function Loadfile_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag','sviewer'),'UserData');
if isempty(Data)
   errordlg('No File is open!', 'Error');
   return;
end
if isfield(Data,'Detection')
    setdefault_detection_stop(Data,1);
end
[file,path]=uigetfile('*.*','Open Detection');
if isequal(file,0) | isequal(path,0)
    return;
else
    tmpflag = isequal(Data.HDR.TYPE,'alpha') & isequal(file,'marker');
    tmpflag = tmpflag | isequal(file,Data.File.file);
    tmpflag = tmpflag & isequal(path,Data.File.path);
    if tmpflag,
        H = Data.HDR;
    else
        H = sopen([path file]);
        H = sclose(H);
    end
    try
        typ = H.EVENT.TYP;
        pos = H.EVENT.POS;
        if isfield(H.EVENT,'CHN')
                chn = H.EVENT.CHN;
        else    
                chn = zeros(size(pos));
        end;
        if isfield(H.EVENT,'DUR')
                dur = H.EVENT.DUR;
        else    
                dur = zeros(size(pos));
        end;
    catch
        errordlg('The selected Event-File is incorrect!', 'File Error');
        return;
    end
        Data.Detection.EventMatrix = [round(pos) typ chn round(dur) [1:length(typ)]'];
        set(findobj('Tag','sviewer'),'UserData',Data);
        Startdetection_Callback;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to display different events
function Radiobutton_Display_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag','sviewer'),'UserData');
try
	firstevent = min(Data.Detection.EventMatrix(1,2));
	if firstevent == 33824
        firstevent = 1056;
	end
	s = Data.Eventcodes_txt;
	temp = s.GroupValue;
    temp2 = unique([temp;firstevent]);
    if size(temp2,1) > size(temp,1)
        pos_firstevent = find(temp2 == firstevent)-1;
    else
        pos_firstevent = find(temp2 == firstevent);
    end
catch
    pos_firstevent = 1;
end
set(findobj('Tag', 'Detection_radiobutton_detection'), 'Value', 0);
Event_string = sprintf('%s |',Data.Eventcodes_txt.GroupDesc{:});
Event_string = Event_string(1:end-2);
set(findobj('Tag','listbox_Event'),'String',Event_string,'Enable','on','BackgroundColor',[0.847 0.922 0.922]);
set(findobj('Tag','listbox_Event'),'Value',pos_firstevent); %1
set(findobj('Tag','listbox_detchannels'),'Enable','off');
set(findobj('Tag','listbox_Event_detail'),'Value',1,'MAX',1500,'BackgroundColor',[0.847 0.922 0.922]);
Listbox_MainClass_Callback;
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to executed different events
function Radiobutton_Detection_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag','sviewer'),'UserData');
set(findobj('Tag', 'Detection_radiobutton_display'), 'Value', 0);
%%%set(findobj('Tag','listbox_Event'),'Value',1);
set(findobj('Tag','listbox_Event'),'BackgroundColor',[1 0.91 0.941]);
set(findobj('Tag','listbox_Event_detail'),'Value',1,'MAX',1,'BackgroundColor',[1 0.91 0.941]);
set(findobj('Tag','listbox_Event_detail'),'String',Data.Detection.Eventdetail_string);
set(findobj('Tag','listbox_detchannels'),'Enable','on');
Data.ShowDetection.MainClass = get(findobj('Tag','listbox_Event'),'Value');
Data.ShowDetection.SubClass = 'all Events';
set(findobj('Tag','sviewer'),'UserData',Data);
h = Data.SPlot;
numb_plot = max(size(h));
redraw_plots(h,numb_plot,Data,Data.Slider.tsec);
figure(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw a new event
function Drawline_Callback(hObject, eventdata, handles)

det_active = get(findobj('Tag', 'Detection_radiobutton_detection'), 'Value');
if det_active == 0
    errordlg('Activate Detection-Radiobutton!', 'Error');
    return
end

mousebutton = get(gcf,'SelectionType'); 
if isequal(mousebutton,'extend')
    return
end
if isequal(mousebutton,'alt')
    return
end
Data = get(findobj('Tag','sviewer'),'UserData');
pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
startchannel = round(startchannel);
sel_channel = startchannel + get(gca, 'UserData') -1;
%channel_name = Data.Channel{sel_channel};
%pos_channel = strmatch(channel_name,Data.allChannel);
pos_channel = sel_channel; 
slider_step = get(findobj('Tag','Slider1'),'SliderStep');
pos_slider = Data.Slider.Pos;
X = pos_slider / slider_step(1) * max(Data.HDR.SampleRate) * Data.NoS / 2;
Y = X + max(Data.HDR.SampleRate) * Data.NoS;
pos1 = get(gca,'CurrentPoint');
pos = pos1(1)+X;
h = Data.SPlot;
numb_channel = size(Data.Channel);
if ~isfield(Data.Detection,'EventMatrix')
    Data.Detection.EventMatrix=[]; 
end 
if isempty(Data.Patch)
    Data.Patch.line = pos;
    Data.Patch.handleline = get(gca, 'UserData')+startchannel-1;
    Data.Patch.selChannel = sel_channel;
    Data.Patch.firstline = line('XData', [pos-X pos-X], ...
                   'YData', [-99999 99999], ...
                   'EraseMode', 'xor', ...
                   'Color', [1 0.2 0.2]);
else
    startline = Data.Patch.line;
    if startline > pos
        if startline <= Y
            delete(Data.Patch.firstline);
        end
        try delete(Data.Patch.firstline); end
        Data.Patch = [];
        set(findobj('Tag','sviewer'),'UserData',Data);
        errordlg('The starting time must lie before the end time!', 'Error');
        return
    end
    if X < Data.Patch.firstline
        h = get(gca,'Children');
        exist_line = find(h == Data.Patch.firstline);
        if exist_line
            delete(Data.Patch.firstline);
        end
    end

    sel_event = get(findobj('Tag','Event'),'Value');
    sel_eventdetail = get(findobj('Tag','Event_detail'),'Value');
    s = Data.Eventcodes_txt;
    select_string = get(findobj('Tag','listbox_Event_detail'),'String');
    select_value = get(findobj('Tag','listbox_Event_detail'),'Value');
    if length(select_string) > 40
        pos_desc = strmatch(select_string(select_value,1:40),s.CodeDesc);
    else
        pos_desc = strmatch(select_string(select_value,:),s.CodeDesc,'exact');
    end
    det_typ = s.CodeIndex(pos_desc);

    color = get_color(det_typ,Data);
    
    pos_channel1 = [min(pos_channel,Data.Patch.selChannel):1:max(pos_channel,Data.Patch.selChannel)];
    numb_showplots = abs(pos_channel-Data.Patch.selChannel)+1;
    act_showplots = pos_channel1 - startchannel + 1; 
    act_showplots = act_showplots(find(act_showplots > 0));
    if numb_showplots == Data.HDR.NS
        pos_channel1 = 0;
    end
    for z = 1:numb_showplots
        i = size(Data.Detection.EventMatrix,1);
        if i == 0
            id = 1;
        else
            id = Data.Detection.EventMatrix(end,5)+1;
        end
        Data.Detection.Update = [];
        if sum(pos_channel1) == 0 & sum(pos_channel1) >= 0
            Data.Detection.Update = [round(startline),det_typ,0,round(pos-startline),id];
            Data.Detection.EventMatrix(i+1,:) = [round(startline),det_typ,0,round(pos-startline),id];
            pos_channel1 = -1;
        elseif sum(pos_channel1) >= 1
            Data.Detection.Update = [round(startline),det_typ,pos_channel1(z),round(pos-startline),id];
            Data.Detection.EventMatrix(i+1,:) = [round(startline),det_typ,pos_channel1(z),round(pos-startline),id];
        end
        if startline < X
            startline1 = 0;
            pos1 = pos - X;
        else
            startline1 = startline - X;
            pos1 = pos - X;
        end
    end
    for z = 1:length(act_showplots)
        Data.actPatch = i;
        Data.Patch.actchannelpos = act_showplots(z);
        patch1 = detpatch('position',[round(startline1) -10000 round(pos1-startline1) 20000], ...
                            'facecolor',color, ...
                            'tag', ['detpatch_',int2str(id)], ...
                            'id', id, ...
                            'data', Data);
        Data.Detection.Patch(:,1) = patch1;
    end
    try delete(Data.Patch.firstline); end
    Data.Patch = [];
end

set(findobj('Tag','sviewer'),'UserData',Data);          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display help page
function Help_Callback(hObject, eventdata, handles)
if strncmpi(computer,'PCWIN',2)
	web(['file:///' which('viewer\help\index.htm')]);
else
	web(['file:///' which('viewer/help/index.htm')]);
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% about
function About_Callback(hObject, eventdata, handles)
helpdlg(sprintf([ 'SViewer \n', ...
      'Version 1.52\n\n' ...
      'Copyright (C) 2004 Franz Einspieler <znarfi5@hotmail.com> and\n', ...
      'Copyright (C) 2004, 2008 Alois Schloegl   <a.schloegl@ieee.org>\n\n' ...
      '        University of Technology Graz, Austria\n\n' ...
      'Comments or suggestions may be sent to the author.\n', ...
      'This Software is subject to the GNU General Public License.']), ...
    'About SViewer');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortcuts
function KeyPress_Callback(hObject, eventdata, handles)
switch get(gcf, 'CurrentKey')
    case 'f1'
        Help_Callback;
    case 'f2'
        OpenFile_Callback;
    case 'f3'
        DisplayMenu_Callback;
% % %     case 'f4'
% % %         ChannelsMenu_Callback;
    case 'f5'
        Startdetection_Callback;
    case 'f6'
        Loadfile_Callback;
    case 'f7'
        FileInfoMenu_Callback;
    case 'f8'
        About_Callback;
    case 'f9'
        Exit_Callback;
    otherwise
         return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%