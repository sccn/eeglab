function varargout = guipas(varargin)
% GUIPAS Application M-file for guifastr.fig
%    FIG = GUIPAS launch guifastr GUI.
%    GUIPAS('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 14-Dec-2004 16:11:40

if nargin == 1  % LAUNCH GUI

	pasguifig = openfig(mfilename,'new');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(pasguifig);
    handles.param.method='';
    handles.param.qrsevent='';
    handles.param.npc=get(handles.npc,'string');
    handles.opstatus=0;
    guidata(pasguifig, handles);
    
    % Intialize pop-up menu
    guipas('qrs_popup_Callback',handles.qrs_popup,[],handles,varargin{1})
    % ---------------------

    % Wait for callbacks to run and window to be dismissed:
	uiwait(pasguifig);
    
	if nargout > 0
		varargout{1} = pasguifig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.

% --- Executes on selection change in qrs_popup.
function varargout = qrs_popup_Callback(h, eventdata, handles, varargin)
if nargin > 3
    set(h,'string',varargin{1});
    handles.param.qrsevent=get(h,'value');
    guidata(h, handles);
else
    handles.param.qrsevent=get(h,'value');
    guidata(h, handles);
end

% --- Executes on button press in rapco_method_check.
function varargout = rapco_method_check_Callback(h, eventdata, handles, varargin)
set(handles.gmean_method_check,'value',0);
set(handles.median_method_check,'value',0);
set(handles.mean_method_check,'value',0);
set(handles.npc_text,'ForeGroundColor',[0 0 0.55]);
set(handles.npc,'Enable','on');
handles.param.method='obs';
handles.param.npc=get(handles.npc,'string');
guidata(h, handles);

% --- Executes on editing npc textbox
function varargout = npc_Callback(h, eventdata, handles, varargin)
handles.param.npc=get(handles.npc,'string');
guidata(h, handles);

% --- Executes on button press in mean_method_check.
function varargout = mean_method_check_Callback(h, eventdata, handles, varargin)
set(handles.gmean_method_check,'value',0);
set(handles.median_method_check,'value',0);
set(handles.rapco_method_check,'value',0);
set(handles.npc_text,'ForeGroundColor',[0.55 0.55 0.55]);
set(handles.npc,'Enable','off');
handles.param.method='mean';
guidata(h, handles);

% --- Executes on button press in gmean_method_check.
function varargout = gmean_method_check_Callback(h, eventdata, handles, varargin)
set(handles.mean_method_check,'value',0);
set(handles.median_method_check,'value',0);
set(handles.rapco_method_check,'value',0);
set(handles.npc_text,'ForeGroundColor',[0.55 0.55 0.55]);
set(handles.npc,'Enable','off');
handles.param.method='gmean';
guidata(h, handles);

% --- Executes on button press in median_method_check.
function varargout = median_method_check_Callback(h, eventdata, handles, varargin)
set(handles.mean_method_check,'value',0);
set(handles.gmean_method_check,'value',0);
set(handles.rapco_method_check,'value',0);
set(handles.npc_text,'ForeGroundColor',[0.55 0.55 0.55]);
set(handles.npc,'Enable','off');
handles.param.method='median';
guidata(h, handles);

% --- Executes on button press in ok_button.
function varargout = ok_button_Callback(h, eventdata, handles, varargin)
handles.param.qrsevent=get(handles.qrs_popup,'value');
if isempty(handles.param.method)
    errordlg('Please select artifact construction method.','fmrib_pas error');
    return;
end
handles.opstatus=1;
guidata(h,handles);
uiresume(handles.pasgui_fig);

% --- Executes on button press in help_button.
function varargout = help_button_Callback(h, eventdata, handles, varargin)
pophelp('pop_fmrib_pas.m');

% --- Executes on button press in cancel_button.
function varargout = cancel_button_Callback(h, eventdata, handles, varargin)
handles.opstatus=0;
guidata(h,handles);
uiresume(handles.pasgui_fig);


% --- Executes during object creation, after setting all properties.
function npc_CreateFcn(hObject, eventdata, handles)
