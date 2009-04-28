function [ varargout ] = eeg_toolbox(command)

% eeg_toolbox - Graphical user interface (GUI) to various EEG/ERP tools
%
% The main gui is the primary store for general parameters 
% and provides access to other tools.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:53 $

% Licence:  GNU GPL, no express or implied warranties
% History:  01/2002, Darren.Weber_at_radiology.ucsf.edu
%           08/2002, Darren.Weber_at_radiology.ucsf.edu
%                    added MRI viewer
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('command','var'),
    command = 'init';
elseif isempty(command),
    command = 'init';
end

command = lower(command);

switch command,
case 'init',
otherwise,
    EEGTOOLBOX = get(gcbf,'Userdata');
end

switch command,
    
case 'init',
    EEGTOOLBOX = init;
    
case 'openerp',
    gui_erp_open(EEGTOOLBOX.p,'init',EEGTOOLBOX.gui);

case 'openerf',
    gui_erf_open(EEGTOOLBOX.p,'init',EEGTOOLBOX.gui);

case 'opencnt',
    gui_cnt_open(EEGTOOLBOX.p,'init',EEGTOOLBOX.gui);
    
case 'openeeg',
    gui_eeg_open(EEGTOOLBOX.p,'init',EEGTOOLBOX.gui);
    
case 'opene',
    gui_elec_open(EEGTOOLBOX.p,'init',EEGTOOLBOX.gui);
    
case 'opent',
    gui_mesh_open(EEGTOOLBOX.p,'init',EEGTOOLBOX.gui);
    
case 'openm',
    mri_toolbox('init');
    
case 'defaultreturn',
    % keep this command so the return action is tidy
    %p = EEGTOOLBOX.p; % done below
    
case 'defaultreset',
    EEGTOOLBOX.p = eeg_toolbox_defaults('create');
    
case 'defaultsave',
    eeg_toolbox_defaults('write',EEGTOOLBOX.p);
    
case 'saveas',
    eeg_toolbox_defaults('write_other',EEGTOOLBOX.p);
    
case 'recent',
    
    % -- get recent files list
    
    recentfiles = eeg_toolbox_recent;
    
    % -- remove current menu items for recent files
    
    if isfield(EEGTOOLBOX.menu,'recentfiles'),
        handleIndex = find(ishandle(EEGTOOLBOX.menu.recentfiles));
        delete(EEGTOOLBOX.menu.recentfiles(handleIndex));
    end
    EEGTOOLBOX.menu.recentfiles = [];
    
    % -- recreate recent files menu items
    
    if and(size(recentfiles,2) == 1, isempty(recentfiles{1})),
        if ishandle(EEGTOOLBOX.menu.recent),
            set(EEGTOOLBOX.menu.recent,'Label','No Recent Files');
        end
    else
        if ishandle(EEGTOOLBOX.menu.recent),
            set(EEGTOOLBOX.menu.recent,'Label','Recent Files');
        end
        
        % -- add recent files to menu and setup their callbacks
        
        for i=1:length(recentfiles),
            if ~isempty(recentfiles{i}),
                EEGTOOLBOX.menu.recentfiles(i) = uimenu(EEGTOOLBOX.menu.recent,...
                    'Label',recentfiles{i},...
                    'Callback',strcat('[recentfiles,p] = eeg_toolbox_recent(''',...
                    recentfiles{i},''',''load''); ',...
                    'EEGTOOLBOX = get(gcbf,''Userdata''); ',...
                    'gui_eeg_open(p,''init'',EEGTOOLBOX.gui); ',...
                    'clear EEGTOOLBOX recentfiles;'));
            end
        end
        
        % -- add recent files clear command
        
        EEGTOOLBOX.menu.recentfiles(i+1) = uimenu(EEGTOOLBOX.menu.recent,...
            'Label','Clear All',...
            'Callback',strcat('eeg_toolbox_recent('''',''clear''); ',...
            'eeg_toolbox(''recent''); '));
    end
    
case 'exit',
    close gcbf;
    
otherwise,
    fprintf('...invalid command to eeg_toolbox\n\n');
    
end


switch command,
case 'exit',
otherwise,
    set(EEGTOOLBOX.gui,'UserData',EEGTOOLBOX);
end

if nargout > 0,
    if isfield(EEGTOOLBOX,'p'),
        if ~isempty(EEGTOOLBOX.p),
            varargout{1} = EEGTOOLBOX.p;
        end
    end
end

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paint the GUI
function [H] = init()

% Parameters are supplied in the file defaultfile.
H.p = eeg_toolbox_defaults('read');

GUIwidth  = 250;
GUIheight = 50;

if exist('ctf_read'),
  name = 'EEG/MEG Toolbox';
else
  name = 'EEG Toolbox';
end

H.gui = figure('Name',name,'Tag','EEG_TOOLBOX',...
    'NumberTitle','off','HandleVisibility','callback',...
    'MenuBar','none');
set(H.gui,'Position',[1 1 GUIwidth GUIheight]);  % Activate GUI Figure
movegui(H.gui, 'center');


license = ['EEG_TOOLBOX, Copyright (C) 2004 Darren L. Weber\n',...
'EEG_TOOLBOX comes with ABSOLUTELY NO WARRANTY; for details\n',...
'see `help gpl`.  This is free software, and you are welcome\n',...
'to redistribute it under certain conditions; see `help gpl`\n',... 
'for details.\n\n'];
fprintf(license);


% -- file menu

H.menu.file_menu = uimenu(H.gui,'Label','File',...
                          'Callback','eeg_toolbox(''recent'');');

H.menu.open_eeg   = uimenu(H.menu.file_menu,'Label','Open EEG');
H.menu.open_erp   = uimenu(H.menu.open_eeg,'Label','Average ERP',...
                           'Callback','p = eeg_toolbox(''openERP'');','Accelerator','v');
H.menu.open_cont  = uimenu(H.menu.open_eeg,'Label','Neuroscan Continuous (CNT)',...
                           'Callback','p = eeg_toolbox(''openCNT'');','Accelerator','c');
H.menu.open_epoch = uimenu(H.menu.open_eeg,'Label','Neuroscan Epochs (EEG)',...
                           'Callback','p = eeg_toolbox(''openEEG'');','Accelerator','e');

if exist('ctf_read'),
  H.menu.open_meg = uimenu(H.menu.file_menu,'Label','Open MEG');
  H.menu.open_erf = uimenu(H.menu.open_meg,'Label','Average ERF',...
                           'Callback','p = eeg_toolbox(''openERF'');','Accelerator','f');
end

H.menu.open_elec = uimenu(H.menu.file_menu,'Label','Open Sensors/Electrodes',...
    'Callback','p = eeg_toolbox(''openE'');','Accelerator','s');
H.menu.open_tess = uimenu(H.menu.file_menu,'Label','Open Tesselation',...
    'Callback','p = eeg_toolbox(''openT'');','Accelerator','t');
if exist('avw_img_read.m') == 2,
    H.menu.open_mri = uimenu(H.menu.file_menu,'Label','Open MRI',...
        'Callback','mri = mri_toolbox(''init'');','Accelerator','m');
end

H.menu.recent = uimenu(H.menu.file_menu,'Label','Recent');
H.menu.quit   = uimenu(H.menu.file_menu,'Label','Exit',...
    'Callback','eeg_toolbox(''exit'');','Accelerator','x');

% -- Parameters menu

H.menu.p_menu = uimenu(H.gui,'Label','Parameters');
H.menu.show   = uimenu(H.menu.p_menu,'Label','Return to Workspace',...
    'Callback','p = eeg_toolbox(''defaultreturn'')');
H.menu.reset  = uimenu(H.menu.p_menu,'Label','Reset to Defaults',...
    'Callback','p = eeg_toolbox(''defaultreset'');');
H.menu.save   = uimenu(H.menu.p_menu,'Label','Save Defaults',...
    'Callback','p = eeg_toolbox(''defaultsave'');');
H.menu.saveas = uimenu(H.menu.p_menu,'Label','Save As Data Workspace',...
    'Callback','p = eeg_toolbox(''saveas'');');

% -- help menu

H.menu.Help = uimenu(H.gui,'Label','Help');
H.menu.help = uimenu(H.menu.Help,'Label','Help','Callback','doc eeg_toolbox;');
H.menu.gpl  = uimenu(H.menu.Help,'Label','GPL',...
    'Callback',['eegpath = eeg_toolbox_path; cd(eegpath);',...
        'fid = fopen([eegpath,''\gpl.txt'']);',...
        'fseek(fid,0,''eof''); eof_bytes = ftell(fid); fseek(fid,0,''bof'');',...
        'while ftell(fid) < eof_bytes,',...
        '      tline = fgetl(fid); disp(tline);',...
        'end;']);

return
