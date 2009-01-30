function varargout = guifastr(varargin)
% GUIFASTR Application M-file for guifastr.fig
%    FIG = GUIFASTR launch guifastr GUI.
%    GUIFASTR('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 14-Dec-2004 12:21:18

if nargin == 1  % LAUNCH GUI

	fastrguifig = openfig(mfilename,'new');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fastrguifig);
    handles.param.volumes='';
    handles.param.slices='';
    handles.param.lpf=get(handles.lpf,'string');
    handles.param.L=get(handles.L,'string');
    handles.param.window=get(handles.window,'string');
    handles.param.exc=get(handles.exc,'string');
    handles.param.trigcorrect_check=get(handles.trigcorrect_check,'value');
    handles.param.pre_frac=get(handles.pre_frac,'string');
    handles.param.npc=get(handles.npc,'string');
    handles.param.options_check=get(handles.options_check,'value');
    handles.param.strig=1;
    handles.param.anc_chk=0;
    handles.opstatus=0;
    guidata(fastrguifig, handles);
    
    % Intialize pop-up menu
    guifastr('etype_Callback',handles.etype,[],handles,varargin{1})
    % ---------------------

    % Wait for callbacks to run and window to be dismissed:
	uiwait(fastrguifig);
    
	if nargout > 0
		varargout{1} = fastrguifig;
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


% --------------------------------------------------------------------
function varargout = lpf_Callback(h, eventdata, handles, varargin)
handles.param.lpf=get(h,'string');
if isempty(handles.param.lpf)
    handles.param.lpf='0';
end
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = L_Callback(h, eventdata, handles, varargin)
handles.param.L=get(h,'string');
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = window_Callback(h, eventdata, handles, varargin)
handles.param.window=get(h,'string');
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = exc_Callback(h, eventdata, handles, varargin)
handles.param.exc=get(h,'string');
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = anc_chk_Callback(h, eventdata, handles, varargin)
handles.param.anc_chk=get(h,'value');
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = vol_trig_chk_Callback(h, eventdata, handles, varargin)
switch get(h,'value')
case 1
    set(handles.slice_trig_chk,'value',0);
case 0
    set(handles.slice_trig_chk,'value',1);
end
handles.param.strig=0;
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = slice_trig_chk_Callback(h, eventdata, handles, varargin)
switch get(h,'value')
case 1
    set(handles.vol_trig_chk,'value',0);
case 0
    set(handles.vol_trig_chk,'value',1);
end
handles.param.strig=1;
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = etype_Callback(h, eventdata, handles, varargin)
if nargin > 3
    set(h,'string',varargin{1});
    handles.param.etype=get(h,'value');
    guidata(h, handles);
else
    handles.param.etype=get(h,'value');
    guidata(h, handles);
end

% --------------------------------------------------------------------
function varargout = trigcorrect_check_Callback(h, eventdata, handles, varargin)
handles.param.trigcorrect_check=get(h,'value');
if handles.param.trigcorrect_check==1
    set(handles.volumes_text,'ForeGroundColor',[0 0 0.55]);
    set(handles.slices_text,'ForeGroundColor',[0 0 0.55]);
    set(handles.volumes,'enable','on');
    set(handles.slices,'enable','on');
elseif handles.param.trigcorrect_check==0
    set(handles.volumes_text,'ForeGroundColor',[0.55 0.55 0.55]);
    set(handles.slices_text,'ForeGroundColor',[0.55 0.55 0.55]);
    set(handles.volumes,'enable','inactive');
    set(handles.slices,'enable','inactive');
    set(handles.volumes,'string','');
    set(handles.slices,'string','');
end
    
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = volumes_Callback(h, eventdata, handles, varargin)
handles.param.volumes=get(h,'string');
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = slices_Callback(h, eventdata, handles, varargin)
handles.param.slices=get(h,'string');
guidata(h, handles);
 
% --------------------------------------------------------------------
function varargout = options_check_Callback(h, eventdata, handles, varargin)
handles.param.options_check=get(h,'value');
if handles.param.options_check==1
    set(handles.pre_frac_text,'visible','on');
    set(handles.pre_frac,'visible','on');
    set(handles.npc_text,'visible','on');
    set(handles.npc_note,'visible','on');
    set(handles.npc,'visible','on');
else
    set(handles.pre_frac_text,'visible','off');
    set(handles.pre_frac,'visible','off');
    set(handles.npc_text,'visible','off');
    set(handles.npc_note,'visible','off');
    set(handles.npc,'visible','off');
end    
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = pre_frac_Callback(h, eventdata, handles, varargin)
handles.param.pre_frac=get(h,'string');
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = npc_Callback(h, eventdata, handles, varargin)
handles.param.npc=get(h,'string');
if isempty(handles.param.npc)
    handles.param.npc='auto';
end
    
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = help_button_Callback(h, eventdata, handles, varargin)
pophelp('pop_fmrib_fastr.m');

% --------------------------------------------------------------------
function varargout = cancel_button_Callback(h, eventdata, handles, varargin)
handles.opstatus=0;
guidata(h,handles);
uiresume(handles.guifastr_fig);

% --------------------------------------------------------------------
function varargout = ok_button_Callback(h, eventdata, handles, varargin)
if handles.param.trigcorrect_check==1
	if isempty(handles.param.volumes)
        errordlg('Please enter the number of FMRI volumes','FASTR error');
        return;
	elseif isempty(handles.param.slices)
        errordlg('Please enter the number of slices per FMRI volume','FASTR error');
        return;
	end
else
    handles.param.slices='0';
    handles.param.volumes='0';
end

fprintf('\n-----------------------------------\n');
fprintf('FASTR Settings:\n');
fprintf('---------------\n');
if str2num(handles.param.lpf)>0   
    fprintf('Low-pass filter cutoff frequency (Hz): %d\n',...
                str2num(handles.param.lpf));    
else
    if handles.param.anc_chk==1
        fprintf('Low-pass filter cutoff frequency (Hz): 70\n');
    else
        fprintf('Low-pass filter cutoff frequency (Hz): NO LPF\n');
    end
end
fprintf('Averaging window length: %d\n',floor(str2num(handles.param.window)/2)*2+1);
evalue=get(handles.etype,'value');
etypes=get(handles.etype,'string');
etype=char(etypes{evalue});
fprintf('Artifact-timing events: %s\n',etype);
if handles.param.strig==1
    fprintf('   Events are SLICE events\n');
else
    fprintf('   Events are VOLUME/SECTION events\n');
end
if isempty(handles.param.exc)
    fprintf('All channels will be processed as EEG channels\n');
    fprintf('   WARNING: FASTR might degrade the quality of any non-EEG channels\n');
else
    fprintf('These non-EEG channels will not be processed with OBS: %s\n',handles.param.exc);
end
if handles.param.trigcorrect_check==1
    fprintf('FASTR will attempt to correct for missing slice triggers for:\n');
    fprintf('   FMRI Volumes: %d\n',str2num(handles.param.volumes));
    fprintf('   Slices per volume: %d\n',str2num(handles.param.slices));
else
    fprintf('FASTR will not correct for missing slice triggers.\n');
end
if handles.param.options_check==1
    fprintf('Advanced Options:\n');
    fprintf('   Trigger location within slice artifact (0-1): %s\n',handles.param.pre_frac);
    fprintf('   Number of PCs for use in OBS: %s\n',handles.param.npc);
end
handles.opstatus=1;
guidata(h,handles);
uiresume(handles.guifastr_fig);

% --------------------------------------------------------------------
function pre_frac_CreateFcn(h, eventdata, handles, varargin)
