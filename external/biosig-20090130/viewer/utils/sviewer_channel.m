function varargout = sviewer_channel(varargin)
% SVIEWER_CHANNEL
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
                   'gui_OpeningFcn', @sviewer_channel_OpeningFcn, ...
                   'gui_OutputFcn',  @sviewer_channel_OutputFcn, ...
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
% --- Executes just before sviewer_channel is made visible.
function sviewer_channel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sviewer_channel (see VARARGIN)

% Choose default command line output for sviewer_channel
handles.output = hObject;
set(gcf,'Color',[0.949,0.949,1]);
% Update handles structure
guidata(hObject, handles);

Data = get(findobj('Tag','sviewer'),'UserData');
if isempty(Data)
    return;
else
    Data = get(findobj('Tag', 'sviewer'), 'UserData');
    set(findobj('Tag', 'listbox_file'), 'String', Data.Channel(:,1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in button_move_up.
function button_move_up_Callback(hObject, eventdata, handles)

Data = get(findobj('Tag','sviewer'),'UserData');
allchannels = get(findobj('Tag','listbox_file'),'String');
act_pos = get(findobj('Tag','listbox_file'),'Value');
disp_min = Data.ChannelConf.Display_min;
disp_max = Data.ChannelConf.Display_max;
if act_pos > 1
    neworder = [allchannels(1:act_pos-2);allchannels(act_pos);allchannels(act_pos-1);allchannels(act_pos+1:end)];
    new_min = [disp_min(1:act_pos-2);disp_min(act_pos);disp_min(act_pos-1);disp_min(act_pos+1:end)];
    new_max = [disp_max(1:act_pos-2);disp_max(act_pos);disp_max(act_pos-1);disp_max(act_pos+1:end)];
    Data.ChannelConf.Display_min = new_min;
    Data.ChannelConf.Display_max = new_max;
    set(findobj('Tag', 'listbox_file'), 'String', neworder);
    set(findobj('Tag','listbox_file'),'Value',act_pos-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in button_move_down.
function button_move_down_Callback(hObject, eventdata, handles)
% hObject    handle to button_move_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Data = get(findobj('Tag','sviewer'),'UserData');
allchannels = get(findobj('Tag','listbox_file'),'String');
act_pos = get(findobj('Tag','listbox_file'),'Value');
disp_min = Data.ChannelConf.Display_min;
disp_max = Data.ChannelConf.Display_max;
if act_pos < size(allchannels,1)
    neworder = [allchannels(1:act_pos-1);allchannels(act_pos+1);allchannels(act_pos);allchannels(act_pos+2:end)];
    new_min = [disp_min(1:act_pos-1);disp_min(act_pos+1);disp_min(act_pos);disp_min(act_pos+2:end)];
    new_max = [disp_max(1:act_pos-1);disp_max(act_pos+1);disp_max(act_pos);disp_max(act_pos+2:end)];
    Data.ChannelConf.Display_min = new_min;
    Data.ChannelConf.Display_max = new_max;
    set(findobj('Tag', 'listbox_file'), 'String', neworder);
    set(findobj('Tag','listbox_file'),'Value',act_pos+1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in button_OK.
function button_OK_Callback(hObject, eventdata, handles)
% hObject    handle to button_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Data = get(findobj('Tag', 'sviewer'), 'UserData');
if isempty(Data)
    close;
    errordlg('SViewer is not open!', 'Error');
    return;
end

channel = get(findobj('Tag','listbox_file'),'String');
Data.Channel(:,1) = channel;
for i = 1:size(channel,1)
    Data.Channel{i,2} = i;
end
set(findobj('Tag', 'Slider1'), 'Value',0);
set(findobj('Tag','sviewer'),'UserData',Data);
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to button_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data = get(findobj('Tag', 'sviewer'), 'UserData');
Data.changeChannel = 'no';
set(findobj('Tag','sviewer'),'UserData',Data);
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs from this function are returned to the command line.
function varargout = sviewer_channel_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
