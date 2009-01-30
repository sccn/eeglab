function varargout = sviewer_display(varargin)
% SVIEWER_DISPLAY
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
                   'gui_OpeningFcn', @sviewer_display_OpeningFcn, ...
                   'gui_OutputFcn',  @sviewer_display_OutputFcn, ...
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
% Executes just before sviewer_display is made visible.
function sviewer_display_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sviewer_display (see VARARGIN)

% Choose default command line output for sviewer_display
handles.output = hObject;
set(gcf,'Color',[0.949,0.949,1]);
guidata(hObject, handles);

% activated radiobutton
set(findobj('Tag', 'radiobutton_display'), 'Value', 1);
set(findobj('Tag', 'gototime'), 'Enable', 'off');

Data=get(findobj('Tag','sviewer'),'UserData');
if isempty(Data)
    return;
else
    set(findobj('Tag', 'act_display'), 'String', Data.NoS);
    set(findobj('Tag', 'totallength'), 'String', Data.Total_length_sec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs from this function are returned to the command line.
function varargout = sviewer_display_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newdisplaytime_Callback(hObject, eventdata, handles)

button_OK_Callback;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gototime_Callback(hObject, eventdata, handles)

button_OK_Callback;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in radiobutton_display.
function radiobutton_display_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Value',1);
set(findobj('Tag', 'radiobutton_goto'), 'Value', 0);
set(findobj('Tag', 'newdisplaytime'), 'Enable', 'on');
set(findobj('Tag', 'gototime'), 'Enable', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in radiobutton_goto.
function radiobutton_goto_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_goto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(findobj('Tag', 'radiobutton_display'), 'Value', 0);
set(hObject, 'Value',1);
set(findobj('Tag', 'newdisplaytime'), 'Enable', 'off');
set(findobj('Tag', 'gototime'), 'Enable', 'on');

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

Data.Display.newdisplaytime = '';
Data.Display.gotime = '';
maxtime = Data.Total_length_sec;
file = Data.File.file;
path = Data.File.path;

displaytime = get(findobj('Tag', 'newdisplaytime'), 'String');
gotime = get(findobj('Tag', 'gototime'), 'String');

display_on = (get(findobj('Tag', 'radiobutton_display'),'Value'));
goto_on = (get(findobj('Tag', 'radiobutton_goto'),'Value'));

if display_on
    if isempty(displaytime)
        errordlg('No value for new Displayed-Time!', 'Error');
        return;
    end
    displaytime = str2double(displaytime);
    if displaytime > maxtime
        errordlg(['The value for Displayed-Time is too large ! [max. ' num2str(maxtime) ' sec.]'], 'Error');
        return;
    end
    if isnan(displaytime)
        errordlg('The value for Displayed-Time is Not a Number!', 'Error');
        return;
    end
    if length(displaytime) > 1 | length(displaytime) == 0
        errordlg('The value for Display-Time is Not a Number!', 'Error');
        return;
    end
    if displaytime == 0
        errordlg('The value must be larger than 0!', 'Error');
        return;
    end
    Data.Display.newdisplaytime = displaytime;
end

if goto_on
    if isempty(gotime)
        errordlg('No value for Go To Second!', 'Error');
        return;
    end
    gotime = str2double(gotime);
    if gotime > maxtime
        errordlg(['The value for Go To Second is too large ! [max. ' int2str(maxtime) ' sec.]'], 'Error');
        return;
    end
    if isnan(gotime)
        errordlg('The value for Go To Second is Not a Number!', 'Error');
        return;
    end
    if gotime >= (maxtime-Data.NoS)
        gotime = gotime - Data.NoS;
    end
    Data.Display.gotime = gotime;
end
set(findobj('Tag','sviewer'),'UserData',Data);
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in button_Cancel.
function button_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to button_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data = get(findobj('Tag', 'sviewer'), 'UserData');
Data.Display = [];
set(findobj('Tag','sviewer'),'UserData',Data);
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
