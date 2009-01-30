function [ctf,GUI] = ctf_read_gui(command);

% ctf_read_gui - GUI interface to read data from a CTF .ds folder
%
% [ctf,FIG] = ctf_read_gui( [command] );
% 
% eg,
%     ctf = ctf_read_gui;
% 
% ctf struct has fields:
%
% ctf.folder
% ctf.header
% ctf.setup
% ctf.sensor
% ctf.data
%
% This function calls:
% ctf_read_res4 - to read in header, gain/offset, and sensor information
% ctf_read_meg4 - to read in the data
% 
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %
%      <                                                      > %  
%      <                    DISCLAIMER:                       > %
%      <                                                      > %
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. > %
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   > %
%      <                    OFFICIAL USE.                     > %
%      <                                                      > %
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<> %
%



% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

% Copyright (C) 2004  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Modified: 11/2003, Darren.Weber_at_radiology.ucsf.edu
%                    - modified from NIH code
%                      simply to allocate data into one large struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('command','var'), command = 'init'; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paint the GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch command,
  
  case 'init',
    
    [ctf,GUI] = INIT;
    
  case 'update',
    
    CTFopen = get(gcbf,'Userdata');
    
    ctf = CTFopen.ctf;
    
    % update the folder
    set(CTFopen.handles.EctfFolder,'String',ctf.folder);
    
    CTFopen.CHANNELS = 1:ctf.setup.number_channels;
    CTFopen.TIME     = ctf.setup.time_sec;
    CTFopen.TRIALS   = 1:ctf.setup.number_trials;
    
    % update the channels
    channelvalue = get(CTFopen.handles.PctfCHANNELS,'Value');
    channeltype = CTFopen.channeltypes{channelvalue};
    switch channeltype,
      case 'all',   CTFopen.CHANNELS = 1:ctf.setup.number_channels;
      case 'eeg',   CTFopen.CHANNELS = ctf.sensor.index.eeg;
      case 'meg',   CTFopen.CHANNELS = ctf.sensor.index.meg;
      case 'ref',   CTFopen.CHANNELS = ctf.sensor.index.ref;
      case 'other', CTFopen.CHANNELS = ctf.sensor.index.other;
      otherwise,    CTFopen.CHANNELS = [];
    end;
    set(CTFopen.handles.EctfCHANNELS,'value',CTFopen.CHANNELS);
    set(CTFopen.handles.EctfCHANNELS,'string',num2str(CTFopen.CHANNELS));
    
    % update time
    start = sprintf('%6.5f', CTFopen.TIME(1));
    step  = sprintf('%6.5f',[CTFopen.TIME(2) - CTFopen.TIME(1)]);
    stop  = sprintf('%6.5f', CTFopen.TIME(end));
    timestring = [start,':',step,':',stop];
    
    set(CTFopen.handles.TctfTIME,'String',timestring);
    
    set(CTFopen.handles.EctfTIME,'Value',CTFopen.TIME);
    set(CTFopen.handles.EctfTIME,'String',num2str(CTFopen.TIME'));
    
    % update trials
    start = num2str(CTFopen.TRIALS(1));
    stop =  num2str(CTFopen.TRIALS(end));
    trialstring = [start,':',stop];
    set(CTFopen.handles.TctfTRIALS,'String',trialstring);
    
    set(CTFopen.handles.EctfTRIALS,'Value',CTFopen.TRIALS);
    
    set(CTFopen.handles.EctfTRIALS,'String',num2str(CTFopen.TRIALS));
    set(CTFopen.handles.EctfTRIALS,'Value',CTFopen.TRIALS);
    
    set(gcbf,'Userdata',CTFopen);
    
    
    
  case 'return',
    
    CTFopen = get(gcbf,'Userdata');
    
    ctf = CTFopen.ctf;
    
    CHANNELS = CTFopen.CHANNELS;
    TIME = CTFopen.TIME;
    TRIALS = CTFopen.TRIALS;
    
    ctf = ctf_read_meg4(ctf.folder, ctf, CHANNELS, TIME, TRIALS );
    
    close gcbf;
    
    
  case 'plot',
    
    fprintf('\nCTF_GUI_READ: Plot not implemented yet.\n');
    
%     CTFopen = get(gcbf,'Userdata');
%     
%     CTFopen.p = ctf_open(CTFopen.p,CTFopen.gui);
%     set(CTFopen.gui,'Userdata',CTFopen); % update GUI with CNT data
%     
%     p = gui_updateparent(CTFopen,0);
%     
%     if isequal(get(CTFopen.handles.Bhold,'Value'),0),
%       close gcbf;
%       if isfield(CTFopen,'parent'),
%         parent = CTFopen.parent.gui;
%       else
%         parent = [];
%       end
%     else
%       parent = CTFopen.gui;
%     end
%     
%     plotfig = figure('Name',CTFopen.p.ctf.file,...
%       'NumberTitle','off',...
%       'UserData',CTFopen.p);
%     movegui(plotfig,'center');
%     
%     plot(CTFopen.p.ctf.timeArray,CTFopen.p.ctf.data); axis tight;
%     eeg_plot_metric;
%     
%     if isempty(parent), [Xpoint,Ypoint] = eeg_crosshair('init',CTFopen.p);
%     else                [Xpoint,Ypoint] = eeg_crosshair('init',CTFopen.p,parent);
%     end
    
  case 'save',
    
    fprintf('\nCTF_GUI_READ: Save not implemented yet.\n');
    
  otherwise,
    
    CTFopen = get(gcbf,'Userdata');
    GUI.parent = CTFopen.parent;
    gui_updateparent(GUI);
    close gcbf;
    
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ctf,GUI] = INIT
% GUI General Parameters

ctf = ctf_folder;
ctf = ctf_read_res4(ctf.folder);

% check if the data is greater than 500 Mb
data_size = ctf.setup.number_samples * ctf.setup.number_channels * ctf.setup.number_trials;
data_bytes = data_size * 8;
if data_bytes > 5e9, warning('data is greater than 500 Mb'); end
clear data_size data_bytes;

CTFopen.CHANNELS = 1:ctf.setup.number_channels;
CTFopen.TIME     = ctf.setup.time_sec;
CTFopen.TRIALS   = 1:ctf.setup.number_trials;

GUIwidth  = 800;
GUIheight = 150;

ver = '$Revision: 1.1 $';
name = sprintf('CTF Open [v %s]\n',ver(11:15));

GUI = figure('Name',name,'Tag','CTF_OPEN',...
  'NumberTitle','off',...
  'MenuBar','none','Position',[1 1 GUIwidth GUIheight]);
movegui(GUI,'center');

%Font.FontName   = 'Helvetica';
Font.FontUnits  = 'Pixels';
Font.FontSize   = 12;
Font.FontWeight = 'normal';
Font.FontAngle  = 'normal';

boxdepth = 0.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CTF Data Selection and Parameters

Font.FontWeight = 'bold';

% BROWSE for a .ds folder
browsecommand = strcat('CTFopen = get(gcbf,''Userdata'');',...
  'ctf = CTFopen.ctf; ',...
  'ctf = ctf_folder([],ctf); ',...
  'ctf = ctf_read_res4(ctf.folder); ',...
  'CTFopen.ctf = ctf; ',...
  'set(CTFopen.handles.EctfFolder,''String'',ctf.folder); ',...
  'set(gcbf,''Userdata'',CTFopen);',...
  'clear CTFopen; ctf = ctf_read_gui(''update'');');
G.BctfFile = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized',Font, ...
  'Position',[.01 .80 .17 boxdepth], 'String','CTF Folder (*.ds)',...
  'BackgroundColor',[0.8 0.8 0.0], 'ForegroundColor', [1 1 1], ...
  'HorizontalAlignment', 'center', ...
  'TooltipString','Browse for CTF folder', ...
  'Callback', browsecommand );

Font.FontWeight = 'normal';

% EDIT the .ds folder name
editcommand = strcat('CTFopen = get(gcbf,''Userdata'');',...
  'CTFopen.ctf.folder = get(CTFopen.handles.EctfFolder,''String''); ',...
  'ctf = CTFopen.ctf; ',...
  'ctf = ctf_folder(ctf.folder,ctf); ',...
  'ctf = ctf_read_res4(ctf.folder); ',...
  'CTFopen.ctf = ctf; ',...
  'set(CTFopen.handles.EctfFolder,''String'',ctf.folder); ',...
  'set(gcbf,''Userdata'',CTFopen);',...
  'clear CTFopen; ctf = ctf_read_gui(''update'');');
G.EctfFolder = uicontrol('Parent',GUI,'Style','edit','Units','Normalized',Font,  ...
  'Position',[.20 .80 .75 boxdepth], 'String',ctf.folder,...
  'Callback', editcommand );



%---------------------------------------------------------------
% Channels
height = 0.6;

G.Title_channels = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
  'Position',[.01 height .17 boxdepth],'HorizontalAlignment', 'center','String','Channels:');

% Channel types
CTFopen.channeltypes = {'all', 'eeg', 'meg', 'ref', 'other', '[array]'};
channelpopup = strcat('CTFopen = get(gcbf,''Userdata'');',...
  'channelvalue = get(CTFopen.handles.PctfCHANNELS,''Value''); ',...
  'channeltype = CTFopen.channeltypes{channelvalue}; ',...
  'switch channeltype, ',...
  '  case ''all'',   CTFopen.CHANNELS = 1:CTFopen.ctf.setup.number_channels; ',...
  '  case ''eeg'',   CTFopen.CHANNELS = CTFopen.ctf.sensor.index.eeg; ',...
  '  case ''meg'',   CTFopen.CHANNELS = CTFopen.ctf.sensor.index.meg; ',...
  '  case ''ref'',   CTFopen.CHANNELS = CTFopen.ctf.sensor.index.ref; ',...
  '  case ''other'', CTFopen.CHANNELS = CTFopen.ctf.sensor.index.other; ',...
  '  otherwise,      CTFopen.CHANNELS = []; ',...
  'end; ',...
  'set(CTFopen.handles.EctfCHANNELS,''Value'',CTFopen.CHANNELS); ',...
  'set(CTFopen.handles.EctfCHANNELS,''string'',num2str(CTFopen.CHANNELS)); ',...
  'set(gcbf,''Userdata'',CTFopen);',...
  'clear CTFopen channelvalue channeltype;');
G.PctfCHANNELS = uicontrol('Parent',GUI,'Style','popup','Units','Normalized',Font, ...
  'Position',[.20 height .20 boxdepth], ...
  'String',CTFopen.channeltypes,...
  'Value',1, ...
  'TooltipString','Channel selection - enter values for [array]', ...
  'Callback', channelpopup );

% Channels array
channelarray = strcat('CTFopen = get(gcbf,''Userdata'');',...
  'CTFopen.CHANNELS = str2num(get(CTFopen.handles.EctfCHANNELS,''String'')); ',...
  'set(CTFopen.handles.EctfCHANNELS,''Value'',CTFopen.CHANNELS); ',...
  'set(gcbf,''Userdata'',CTFopen); clear CTFopen;');
G.EctfCHANNELS = uicontrol('Parent',GUI,'Style','edit','Units','Normalized',Font, ...
  'Position',[.40 height .50 boxdepth], ...
  'String',num2str(CTFopen.CHANNELS),'Value',CTFopen.CHANNELS,...
  'TooltipString','Channels array  - enter values for [array] selection', ...
  'Callback', channelarray );
G.EctfCHANNELSCLEAR = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized',Font, ...
  'Position',[.90 height .05 boxdepth], ...
  'String','Clear',...
  'TooltipString','Clear channels array', ...
  'Callback', strcat('CTFopen = get(gcbf,''Userdata'');',...
  'CTFopen.CHANNELS = []; ',...
  'set(CTFopen.handles.EctfCHANNELS,''String'',''''); ',...
  'set(CTFopen.handles.EctfCHANNELS,''Value'',CTFopen.CHANNELS); ',...
  'set(gcbf,''Userdata'',CTFopen); clear CTFopen;') );




%---------------------------------------------------------------
% Time
height = 0.4;

G.Title_channels = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
  'Position',[.01 height .17 boxdepth],'HorizontalAlignment', 'center','String','Time (sec):');

% Time array
start = sprintf('%6.5f', CTFopen.TIME(1));
step  = sprintf('%6.5f',[CTFopen.TIME(2) - CTFopen.TIME(1)]);
stop  = sprintf('%6.5f', CTFopen.TIME(end));
timestring = [start,':',step,':',stop];

G.TctfTIME = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
  'Position',[.20 height .20 boxdepth],'HorizontalAlignment', 'center',...
  'TooltipString','Time array (sec)', ...
  'String',timestring);


timearray = strcat('CTFopen = get(gcbf,''Userdata'');',...
  'CTFopen.TIME = str2num(get(CTFopen.handles.EctfTIME,''String''))''; ',...
  'set(CTFopen.handles.EctfTIME,''Value'',CTFopen.TIME); ',...
  'set(gcbf,''Userdata'',CTFopen); clear CTFopen;');
G.EctfTIME = uicontrol('Parent',GUI,'Style','edit','Units','Normalized',Font, ...
  'Position',[.40 height .50 boxdepth], ...
  'String',num2str(CTFopen.TIME'),'Value',CTFopen.TIME,...
  'TooltipString','Time array (sec)', ...
  'Callback', timearray );
G.EctfTIMECLEAR = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized',Font, ...
  'Position',[.90 height .05 boxdepth], ...
  'String','Clear',...
  'TooltipString','Clear time array', ...
  'Callback', strcat('CTFopen = get(gcbf,''Userdata'');',...
  'CTFopen.TIME = []; ',...
  'set(CTFopen.handles.EctfTIME,''String'',''''); ',...
  'set(CTFopen.handles.EctfTIME,''Value'',CTFopen.TIME); ',...
  'set(gcbf,''Userdata'',CTFopen); clear CTFopen;') );



%---------------------------------------------------------------
% Trials
height = 0.2;

G.Title_trials = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
  'Position',[.01 height .17 boxdepth],'HorizontalAlignment', 'center','String','Trials:');

start = num2str(CTFopen.TRIALS(1));
stop =  num2str(CTFopen.TRIALS(end));
trialstring = [start,':',stop];

G.TctfTRIALS = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
  'Position',[.20 height .20 boxdepth],'HorizontalAlignment', 'center',...
  'String',trialstring);

% Trials array
trialsarray = strcat('CTFopen = get(gcbf,''Userdata'');',...
  'CTFopen.TRIALS = str2num(get(CTFopen.handles.EctfTRIALS,''String'')); ',...
  'set(CTFopen.handles.EctfTRIALS,''Value'',CTFopen.TRIALS); ',...
  'set(gcbf,''Userdata'',CTFopen); clear CTFopen;');
G.EctfTRIALS = uicontrol('Parent',GUI,'Style','edit','Units','Normalized',Font, ...
  'Position',[.40 height .50 boxdepth], ...
  'String',num2str(CTFopen.TRIALS),'Value',CTFopen.TRIALS,...
  'TooltipString','Trials array', ...
  'Callback', trialsarray );
G.EctfTRIALSCLEAR = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized',Font, ...
  'Position',[.90 height .05 boxdepth], ...
  'String','Clear',...
  'TooltipString','Clear trials array', ...
  'Callback', strcat('CTFopen = get(gcbf,''Userdata'');',...
  'CTFopen.TRIALS = []; ',...
  'set(CTFopen.handles.EctfTRIALS,''String'',''''); ',...
  'set(CTFopen.handles.EctfTRIALS,''Value'',CTFopen.TRIALS); ',...
  'set(gcbf,''Userdata'',CTFopen); clear CTFopen;') );


% % PLOT: Load & plot the data!
% G.Bplot = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
%   'Position',[.20 .01 .18 .2],...
%   'String','PLOT','BusyAction','queue',...
%   'TooltipString','Plot the EEG data and return p struct.',...
%   'BackgroundColor',[0.0 0.5 0.0],...
%   'ForegroundColor', [1 1 1], 'HorizontalAlignment', 'center',...
%   'Callback',strcat('CTFopen = get(gcbf,''Userdata'');',...
%   'p = ctf_read_gui(CTFopen.p,''plot'');',...
%   'clear CTFopen;'));
% 
% % Save As
% G.Bsave = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
%   'Position',[.40 .01 .18 .2],'HorizontalAlignment', 'center',...
%   'String','SAVE AS','TooltipString','EEG File Conversion Tool (not implemented yet)',...
%   'BusyAction','queue',...
%   'Visible','off',...
%   'BackgroundColor',[0.0 0.0 0.75],...
%   'ForegroundColor', [1 1 1], 'HorizontalAlignment', 'center',...
%   'Callback',strcat('CTFopen = get(gcbf,''Userdata'');',...
%   'p = ctf_read_gui(CTFopen.p,''save'');',...
%   'clear CTFopen;'));

Font.FontWeight = 'bold';

% Quit, return file parameters
G.Breturn = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
  'Position',[.70 .01 .10 boxdepth],...
  'String','RETURN','BusyAction','queue',...
  'TooltipString','Return p struct to workspace and parent GUI.',...
  'BackgroundColor',[0.0 0.75 0.0],...
  'ForegroundColor', [1 1 1], 'HorizontalAlignment', 'center',...
  'Callback','ctf = ctf_read_gui(''return'');' );

% Cancel
G.Bcancel = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
  'Position',[.85 .01 .10 boxdepth],...
  'String','CANCEL','BusyAction','queue',...
  'TooltipString','Close, do not read data.',...
  'BackgroundColor',[0.75 0.0 0.0],...
  'ForegroundColor', [1 1 1], 'HorizontalAlignment', 'center',...
  'Callback','close(gcbf);');

% Store userdata
CTFopen.gui = GUI;          
CTFopen.handles = G;
CTFopen.ctf = ctf;
set(GUI,'Userdata',CTFopen);
set(GUI,'HandleVisibility','callback');

return
