function varargout = sviewer_channel_conf(varargin)
% SVIEWER_CHANNEL_CONF
% Select HELP in the Info menu 
%
% Version 1.0, November 2004
% Copyright by (C) Franz Einspieler <znarfi5@hotmail.com> and
%                  Alois Schloegl   <a.schloegl@ieee.org>
% University of Technology Graz, Austria
%
% This is part of the BIOSIG-toolbox http://biosig.sf.net/
% Comments or suggestions may be sent to the author.
% This Software is subject to the GNU public license.

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sviewer_channel_conf_OpeningFcn, ...
                   'gui_OutputFcn',  @sviewer_channel_conf_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sviewer_channel_conf_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for sviewer_channel_conf
handles.output = hObject;
set(gcf,'Color',[0.949,0.949,1]);
guidata(hObject, handles);

Data=get(findobj('Tag','sviewer'),'UserData');
if isempty(Data)
    return;
else
    button = Data.Channelconf.whichbutton;
    
    pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
    startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
    startchannel = round(startchannel);
    button = button + startchannel - 1;
    
    channel_name = Data.Channel{button};
    pos = strmatch(channel_name,Data.allChannel);
    set(findobj('Tag', 'text_channel'), 'String',channel_name);
    display_min = Data.ChannelConf.Display_min(pos,:);
    set(findobj('Tag', 'edit_DisplayMin'), 'String',display_min);
    display_max = Data.ChannelConf.Display_max(pos,:);
    set(findobj('Tag', 'edit_DisplayMax'), 'String',display_max);
    try
        dimension = Data.HDR.PhysDim(pos,:);
    catch
        dimension = '[?]';
    end
    set(findobj('Tag', 'text_ph_min'), 'String',dimension);
    set(findobj('Tag', 'text_ph_max'), 'String',dimension);
    set(findobj('Tag', 'text_d_min'), 'String',dimension);
    set(findobj('Tag', 'text_d_max'), 'String',dimension);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = sviewer_channel_conf_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton_OK_Callback(hObject, eventdata, handles)

value_disp = get(findobj('Tag', 'Applyall_disp'), 'Value');

if value_disp == 1
    select = 2;
else
    select = 0;
end
close_figure(select);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close window
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate new range
function pushbutton_autoscale_Callback(hObject, eventdata, handles)

Data=get(findobj('Tag','sviewer'),'UserData');
button = Data.Channelconf.whichbutton;
pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
startchannel = round(startchannel);
button = button + startchannel - 1;
channel_name = Data.Channel{button};
pos = strmatch(channel_name,Data.allChannel);
sample_min = min(Data.signal(:,pos));
sample_max = max(Data.signal(:,pos));
if sample_min == sample_max
    sample_min = sample_min -10;
    sample_max = sample_max +10;
end
set(findobj('Tag', 'edit_DisplayMin'), 'String',sample_min);
set(findobj('Tag', 'edit_DisplayMax'), 'String',sample_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change the display minimum
function edit_DisplayMin_Callback(hObject, eventdata, handles)

% % % Data = get(findobj('Tag', 'sviewer'), 'UserData');
% % % button = Data.Channelconf.whichbutton;
% % % pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
% % % startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
% % % startchannel = round(startchannel);
% % % button = button + startchannel - 1;
% % % display_min = Data.ChannelConf.Display_min(button,:);
% % % display_min = round(display_min*100)/100;
% % % set(findobj('Tag', 'edit_DisplayMin'), 'String',display_min);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%change the display maximum
%%%function edit_DisplayMax_Callback(hObject, eventdata, handles)

% % % Data = get(findobj('Tag', 'sviewer'), 'UserData');
% % % button = Data.Channelconf.whichbutton;
% % % pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
% % % startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
% % % startchannel = round(startchannel);
% % % button = button + startchannel - 1;
% % % display_max = Data.ChannelConf.Display_max(button,:);
% % % display_max = round(display_max*100)/100;
% % % set(findobj('Tag', 'edit_DisplayMax'), 'String',display_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % change the zoom factor
% % % function edit_Zoomfactor_Callback(hObject, eventdata, handles)
% % % 
% % % Data = get(findobj('Tag', 'sviewer'), 'UserData');
% % % scale_factor = Data.ChannelConf.Scale;
% % % set(findobj('Tag', 'edit_Zoomfactor'), 'String',scale_factor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close window
function close_figure(select)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
display_min_t = get(findobj('Tag', 'edit_DisplayMin'), 'String');
display_max_t = get(findobj('Tag', 'edit_DisplayMax'), 'String');
scale_factor_t = get(findobj('Tag', 'zoom_factor'), 'String');
display_min = str2num(display_min_t);
display_max = str2num(display_max_t);
scale_factor = str2num(scale_factor_t);

if length(display_min) > 1 | length(display_min) == 0
    errordlg('The value for Display-Min is Not a Number!', 'Error');
    return;
end
if length(display_max) > 1 | length(display_max) == 0
    errordlg('The value for Display-Max is Not a Number!', 'Error');
    return;
end
if display_min > display_max
    errordlg('Display-Min is larger than Display-Max !', 'Error');
    return;
end

if select == 0
    button = Data.Channelconf.whichbutton;
    pos_sliderchannel = get(findobj('Tag','Slider_Channel'), 'Value');
    startchannel = (size(Data.Channel,1) - Data.NS) - pos_sliderchannel * (size(Data.Channel,1) - Data.NS) + 1;
    startchannel = round(startchannel);
    button = button + startchannel - 1;
    channel_name = Data.Channel{button};
    pos = strmatch(channel_name,Data.allChannel);
    Data.ChannelConf.Display_min(pos,:) = display_min;
    Data.ChannelConf.Display_max(pos,:) = display_max;
    Data.ChannelConf.Scale = scale_factor; 
    set(findobj('Tag','sviewer'),'UserData',Data);
    close;
    return;
end
if select == 2
    Data.ChannelConf.Display_min(1:end,:) = display_min;
    Data.ChannelConf.Display_max(1:end,:) = display_max;
    Data.ChannelConf.Scale = scale_factor; 
    set(findobj('Tag','sviewer'),'UserData',Data);
    close;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in Applyall_disp.
function Applyall_disp_Callback(hObject, eventdata, handles)

set(findobj('Tag', 'Applyall_auto'), 'Value', 0);
set(findobj('Tag', 'edit_DisplayMin'), 'Enable', 'on');
set(findobj('Tag', 'edit_DisplayMax'), 'Enable', 'on');
set(findobj('Tag', 'pushbutton_autoscale'), 'Enable', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




