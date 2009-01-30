function varargout = sviewer_fileinfo(varargin)
% SVIEWER_FILEINFO
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
                   'gui_OpeningFcn', @sviewer_fileinfo_OpeningFcn, ...
                   'gui_OutputFcn',  @sviewer_fileinfo_OutputFcn, ...
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


% Executes just before sviewer_fileinfo is made visible.
function sviewer_fileinfo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sviewer_fileinfo (see VARARGIN)

% Choose default command line output for sviewer_fileinfo
handles.output = hObject;
set(gcf,'Color',[0.949,0.949,1]);
% Update handles structure
guidata(hObject, handles);
read_info(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs from this function are returned to the command line.
function varargout = sviewer_fileinfo_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read file info
function read_info(Data)

Data=Data{1};
i = 1;
if isfield(Data.HDR,'FILE')
    set(findobj('Tag',['text' num2str(i)]),'String','Filename:');
    set(findobj('Tag',['text' num2str(i) '_1']),'String',[Data.HDR.FILE.Name '.' Data.HDR.FILE.Ext]);
    i=i+1;
end
if isfield(Data.HDR,'TYPE')
    set(findobj('Tag',['text' num2str(i)]),'String','Type:');
    set(findobj('Tag',['text' num2str(i) '_1']),'String',Data.HDR.TYPE);
    i=i+1;
end
if isfield(Data.HDR,'SampleRate')
    set(findobj('Tag',['text' num2str(i)]),'String','SampleRate:');
    set(findobj('Tag',['text' num2str(i) '_1']),'String',Data.HDR.SampleRate);
    i=i+1;
end
if isfield(Data.HDR,'NS')
    set(findobj('Tag',['text' num2str(i)]),'String','Nummber of Channels:');
    set(findobj('Tag',['text' num2str(i) '_1']),'String',Data.HDR.NS);
    i=i+1;
end
if isfield(Data.HDR,'PID')
    set(findobj('Tag',['text' num2str(i)]),'String','PID:');
    set(findobj('Tag',['text' num2str(i) '_1']),'String',Data.HDR.PID);
    i=i+1;
end
if isfield(Data.HDR,'RID')
    set(findobj('Tag',['text' num2str(i)]),'String','RID:');
    set(findobj('Tag',['text' num2str(i) '_1']),'String',Data.HDR.RID); 
    i=i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close window
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%