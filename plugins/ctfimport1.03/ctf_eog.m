function [ctf,GUI] = ctf_eog(command);

% ctf_eog - GUI interface to check EOG data in a CTF .ds folder
%
% [ctf,FIG] = ctf_eog( [command] );
% 
% eg,
%     ctf = ctf_eog;
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


% $Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $

% Copyright (C) 2006  Darren L. Weber
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

% Created: 01/2006, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist('command','var'), command = 'init'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paint the GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch command,
  
 case 'init',
  
  [ctf,GUI] = INIT;
  
 case 'update',
  
  CTFeog = get(gcbf,'Userdata');
  
  ctf = CTFeog.ctf;
  
  CTFeog.CHANNELS = 1:ctf.setup.number_channels;
  CTFeog.TIME     = ctf.setup.time_sec;
  CTFeog.TRIALS   = 1:ctf.setup.number_trials;
  
  set(gcbf,'Userdata',CTFeog);
  
  
 case 'return',
  
  CTFeog = get(gcbf,'Userdata');
  
  ctf = CTFeog.ctf;
  
  close gcbf;    
  
 case 'plot',
  
  fprintf('\nCTF_EOG: Plot not implemented yet.\n');
  
 case 'save',
  
  fprintf('\nCTF_EOG: Save not implemented yet.\n');
  
 otherwise,
  
  CTFeog = get(gcbf,'Userdata');
  GUI.parent = CTFeog.parent;
  gui_updateparent(GUI);
  close gcbf;
  
end

return






% -------------------------------------------

fig = figure;

for e = allTrials,
  
  set(fig,'Name',sprintf('Trial %04d',e))
  
  veogData = ctf.data(:,veog,e);
  heogData = ctf.data(:,heog,e);

  veogBaseline = mean(veogData(1:blPoints));
  heogBaseline = mean(heogData(1:blPoints));

  veogDataBaselined = veogData - veogBaseline;
  heogDataBaselined = heogData - heogBaseline;

  %plot(timeSec,[veogDataBaselined,heogDataBaselined]);
  %legend('VEOG','HEOG')

  % EEG differential must be normalized by the sample rate (sec)
  veogDataDiff = [0; diff(veogDataBaselined)] / ctf.setup.sample_sec;
  heogDataDiff = [0; diff(heogDataBaselined)] / ctf.setup.sample_sec;
  
  % check if this is an artifact trial
  EYEBLINKaTrial = find(EYEBLINKaTrials == e);
  EYEBLINKdTrial = find(EYEBLINKdTrials == e);
  EYEMVMTaTrial  = find(EYEMVMTaTrials == e);
  EYEMVMTdTrial  = find(EYEMVMTdTrials == e);
  
  % -------------------------------
  % VEOG check

  
  
  % amplitude threshold
  t = thresholdDetect.veog.amp; % 80 uV
                                % detect values above threshold
  eyeblinkA = find(abs(veogDataBaselined) > t);
  eyeblinkWave = veogDataBaselined * 0;
  eyeblinkWave(eyeblinkA) = max(abs(veogDataBaselined));

  veogAxisAmp = subplot(2,2,1);
  plot(timeSec, [veogDataBaselined, eyeblinkWave]); hold on
  legend('VEOG', 'VEOG Blink (Amp)')
  plot(timeSec, [veogData*0, (veogData*0)+t, (veogData*0)-t],'k:'); hold off
  ylabel('Volt')
  xlabel('Time (sec)')
  set(veogAxisAmp,'XLim',[ctf.setup.start_sec, ctf.setup.end_sec]);
  ampLimit = thresholdDetect.veog.amp + (thresholdDetect.veog.amp * 0.5);
  set(veogAxisAmp,'YLim',[-ampLimit, ampLimit])
  set(veogAxisAmp,'YTick',-ampLimit:10e-6:ampLimit)
  if EYEBLINKaTrial,
    set(veogAxisAmp,'color',[1 0.75 0.75])
  else
    set(veogAxisAmp,'color',[1 1 1])
  end
  
  
  
  % detect derivatives above 20 mV/sec
  t = thresholdDetect.veog.der;
  eyeblinkD = find(abs(veogDataDiff) > t);
  eyeblinkWave = veogDataDiff * 0;
  eyeblinkWave(eyeblinkD) = max(abs(veogDataDiff));
  
  veogAxisDer = subplot(2,2,3);
  plot(timeSec,[veogDataDiff, eyeblinkWave]); hold on
  legend('diff(VEOG)/sec', 'VEOG Blink (Der)')
  plot(timeSec, [veogData*0, (veogData*0)+t, (veogData*0)-t],'k:'); hold off
  ylabel('Volt / sec')
  xlabel('Time (sec)')
  set(veogAxisDer,'XLim',[ctf.setup.start_sec, ctf.setup.end_sec]);
  ampLimit = thresholdDetect.veog.der + (thresholdDetect.veog.der * 0.5);
  set(veogAxisDer,'YLim',[-ampLimit, ampLimit])
  %set(veogAxisDer,'YTick',-ampLimit:10e-6:ampLimit)
  if EYEBLINKdTrial,
    set(veogAxisDer,'color',[1 0.75 0.75])
  else
    set(veogAxisDer,'color',[1 1 1])
  end
  
  % -------------------------------
  % HEOG check

  
  
  % amplitude threshold
  t = thresholdDetect.heog.amp;
  % detect values above threshold
  eyeMoveA = find(abs(heogDataBaselined) > t);
  eyeMoveWave = heogDataBaselined * 0;
  eyeMoveWave(eyeMoveA) = max(abs(heogDataBaselined));
  
  heogAxisAmp = subplot(2,2,2);
  plot(timeSec, [heogDataBaselined, eyeMoveWave]); hold on
  legend('HEOG', 'HEOG Move (Amp)')
  plot(timeSec, [heogData*0, (heogData*0)+t, (heogData*0)-t],'k:'); hold off
  ylabel('Volt')
  xlabel('Time (sec)')
  set(heogAxisAmp,'XLim',[ctf.setup.start_sec, ctf.setup.end_sec]);
  ampLimit = thresholdDetect.heog.amp + (thresholdDetect.heog.amp * 0.5);
  set(heogAxisAmp,'YLim',[-ampLimit, ampLimit])
  set(heogAxisAmp,'YTick',-ampLimit:10e-6:ampLimit)
  if EYEMVMTaTrial,
    set(heogAxisAmp,'color',[1 0.75 0.75])
  else
    set(heogAxisAmp,'color',[1 1 1])
  end


  
  % HEOG movements may accelerate at about 20 mV/sec
  t = thresholdDetect.heog.der;
  eyeMoveD = find(abs(heogDataDiff) > t);
  eyeMoveWave = heogDataDiff * 0;
  eyeMoveWave(eyeMoveD) = max(abs(heogDataDiff));

  heogAxisDer = subplot(2,2,4);
  plot(timeSec, [heogDataDiff, eyeMoveWave]); hold on
  legend('diff(HEOG)/sec', 'HEOG Move (Der)')
  plot(timeSec, [heogData*0, (heogData*0)+t, (heogData*0)-t],'k:'); hold off
  ylabel('Volt / sec')
  xlabel('Time (sec)')
  set(heogAxisDer,'XLim',[ctf.setup.start_sec, ctf.setup.end_sec]);
  ampLimit = thresholdDetect.heog.der + (thresholdDetect.heog.der * 0.5);
  set(heogAxisDer,'YLim',[-ampLimit, ampLimit])
  %set(heogAxisDer,'YTick',-ampLimit:10e-6:ampLimit)
  if EYEMVMTdTrial,
    set(heogAxisDer,'color',[1 0.75 0.75])
  else
    set(heogAxisDer,'color',[1 1 1])
  end
  
  pause(10)

end

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ctf,GUI] = INIT
% GUI General Parameters

  ctf = ctf_folder;
  ctf = ctf_read_res4(ctf.folder);
  
  
  thresholdDetect.veog.amp = 50e-6; % 50 uV
  thresholdDetect.veog.der = 25e-3; % 25 mV/sec
  thresholdDetect.heog.amp = 25e-6; % 25 uV
  thresholdDetect.heog.der = 25e-3; % 25 mV/sec
  
  

  CTFeog.CHANNELS = 1:ctf.setup.number_channels;
  CTFeog.TIME     = ctf.setup.time_sec;
  CTFeog.TRIALS   = 1:ctf.setup.number_trials;

  GUIwidth  = 800;
  GUIheight = 600;

  ver = '$Revision: 1.1 $';
  name = sprintf('CTF EOG [v %s]\n',ver(11:15));

  GUI = figure('Name',name,'Tag','CTF_EOG',...
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
  browsecommand = strcat('CTFeog = get(gcbf,''Userdata'');',...
                         'ctf = CTFeog.ctf; ',...
                         'ctf = ctf_folder([],ctf); ',...
                         'ctf = ctf_read_res4(ctf.folder); ',...
                         'CTFeog.ctf = ctf; ',...
                         'set(CTFeog.handles.EctfFolder,''String'',ctf.folder); ',...
                         'set(gcbf,''Userdata'',CTFeog);',...
                         'clear CTFeog; ctf = ctf_read_gui(''update'');');
  G.BctfFile = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized',Font, ...
                         'Position',[.01 .80 .17 boxdepth], 'String','CTF Folder (*.ds)',...
                         'BackgroundColor',[0.8 0.8 0.0], 'ForegroundColor', [1 1 1], ...
                         'HorizontalAlignment', 'center', ...
                         'TooltipString','Browse for CTF folder', ...
                         'Callback', browsecommand );

  Font.FontWeight = 'normal';

  % EDIT the .ds folder name
  editcommand = strcat('CTFeog = get(gcbf,''Userdata'');',...
                       'CTFeog.ctf.folder = get(CTFeog.handles.EctfFolder,''String''); ',...
                       'ctf = CTFeog.ctf; ',...
                       'ctf = ctf_folder(ctf.folder,ctf); ',...
                       'ctf = ctf_read_res4(ctf.folder); ',...
                       'CTFeog.ctf = ctf; ',...
                       'set(CTFeog.handles.EctfFolder,''String'',ctf.folder); ',...
                       'set(gcbf,''Userdata'',CTFeog);',...
                       'clear CTFeog; ctf = ctf_read_gui(''update'');');
  G.EctfFolder = uicontrol('Parent',GUI,'Style','edit','Units','Normalized',Font,  ...
                           'Position',[.20 .80 .75 boxdepth], 'String',ctf.folder,...
                           'Callback', editcommand );



  %---------------------------------------------------------------
  % Channels
  height = 0.6;

  G.Title_channels = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
                               'Position',[.01 height .17 boxdepth],'HorizontalAlignment', 'center','String','Channels:');

  % Channel types
  CTFeog.channeltypes = {'all', 'eeg', 'meg', 'ref', 'other', '[array]'};
  channelpopup = strcat('CTFeog = get(gcbf,''Userdata'');',...
                        'channelvalue = get(CTFeog.handles.PctfCHANNELS,''Value''); ',...
                        'channeltype = CTFeog.channeltypes{channelvalue}; ',...
                        'switch channeltype, ',...
                        '  case ''all'',   CTFeog.CHANNELS = 1:CTFeog.ctf.setup.number_channels; ',...
                        '  case ''eeg'',   CTFeog.CHANNELS = CTFeog.ctf.sensor.index.eeg; ',...
                        '  case ''meg'',   CTFeog.CHANNELS = CTFeog.ctf.sensor.index.meg; ',...
                        '  case ''ref'',   CTFeog.CHANNELS = CTFeog.ctf.sensor.index.ref; ',...
                        '  case ''other'', CTFeog.CHANNELS = CTFeog.ctf.sensor.index.other; ',...
                        '  otherwise,      CTFeog.CHANNELS = []; ',...
                        'end; ',...
                        'set(CTFeog.handles.EctfCHANNELS,''Value'',CTFeog.CHANNELS); ',...
                        'set(CTFeog.handles.EctfCHANNELS,''string'',num2str(CTFeog.CHANNELS)); ',...
                        'set(gcbf,''Userdata'',CTFeog);',...
                        'clear CTFeog channelvalue channeltype;');
  G.PctfCHANNELS = uicontrol('Parent',GUI,'Style','popup','Units','Normalized',Font, ...
                             'Position',[.20 height .20 boxdepth], ...
                             'String',CTFeog.channeltypes,...
                             'Value',1, ...
                             'TooltipString','Channel selection - enter values for [array]', ...
                             'Callback', channelpopup );

  % Channels array
  channelarray = strcat('CTFeog = get(gcbf,''Userdata'');',...
                        'CTFeog.CHANNELS = str2num(get(CTFeog.handles.EctfCHANNELS,''String'')); ',...
                        'set(CTFeog.handles.EctfCHANNELS,''Value'',CTFeog.CHANNELS); ',...
                        'set(gcbf,''Userdata'',CTFeog); clear CTFeog;');
  G.EctfCHANNELS = uicontrol('Parent',GUI,'Style','edit','Units','Normalized',Font, ...
                             'Position',[.40 height .50 boxdepth], ...
                             'String',num2str(CTFeog.CHANNELS),'Value',CTFeog.CHANNELS,...
                             'TooltipString','Channels array  - enter values for [array] selection', ...
                             'Callback', channelarray );
  G.EctfCHANNELSCLEAR = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized',Font, ...
                                  'Position',[.90 height .05 boxdepth], ...
                                  'String','Clear',...
                                  'TooltipString','Clear channels array', ...
                                  'Callback', strcat('CTFeog = get(gcbf,''Userdata'');',...
                                                    'CTFeog.CHANNELS = []; ',...
                                                    'set(CTFeog.handles.EctfCHANNELS,''String'',''''); ',...
                                                    'set(CTFeog.handles.EctfCHANNELS,''Value'',CTFeog.CHANNELS); ',...
                                                    'set(gcbf,''Userdata'',CTFeog); clear CTFeog;') );




  %---------------------------------------------------------------
  % Time
  height = 0.4;

  G.Title_channels = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
                               'Position',[.01 height .17 boxdepth],'HorizontalAlignment', 'center','String','Time (sec):');

  % Time array
  start = sprintf('%6.5f', CTFeog.TIME(1));
  step  = sprintf('%6.5f',[CTFeog.TIME(2) - CTFeog.TIME(1)]);
  stop  = sprintf('%6.5f', CTFeog.TIME(end));
  timestring = [start,':',step,':',stop];

  G.TctfTIME = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
                         'Position',[.20 height .20 boxdepth],'HorizontalAlignment', 'center',...
                         'TooltipString','Time array (sec)', ...
                         'String',timestring);


  timearray = strcat('CTFeog = get(gcbf,''Userdata'');',...
                     'CTFeog.TIME = str2num(get(CTFeog.handles.EctfTIME,''String''))''; ',...
                     'set(CTFeog.handles.EctfTIME,''Value'',CTFeog.TIME); ',...
                     'set(gcbf,''Userdata'',CTFeog); clear CTFeog;');
  G.EctfTIME = uicontrol('Parent',GUI,'Style','edit','Units','Normalized',Font, ...
                         'Position',[.40 height .50 boxdepth], ...
  'String',num2str(CTFeog.TIME'),'Value',CTFeog.TIME,...
  'TooltipString','Time array (sec)', ...
  'Callback', timearray );
G.EctfTIMECLEAR = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized',Font, ...
  'Position',[.90 height .05 boxdepth], ...
  'String','Clear',...
  'TooltipString','Clear time array', ...
  'Callback', strcat('CTFeog = get(gcbf,''Userdata'');',...
  'CTFeog.TIME = []; ',...
  'set(CTFeog.handles.EctfTIME,''String'',''''); ',...
  'set(CTFeog.handles.EctfTIME,''Value'',CTFeog.TIME); ',...
  'set(gcbf,''Userdata'',CTFeog); clear CTFeog;') );



%---------------------------------------------------------------
% Trials
height = 0.2;

G.Title_trials = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
  'Position',[.01 height .17 boxdepth],'HorizontalAlignment', 'center','String','Trials:');

start = num2str(CTFeog.TRIALS(1));
stop =  num2str(CTFeog.TRIALS(end));
trialstring = [start,':',stop];

G.TctfTRIALS = uicontrol('Parent',GUI,'Style','text','Units','Normalized',Font, ...
  'Position',[.20 height .20 boxdepth],'HorizontalAlignment', 'center',...
  'String',trialstring);

% Trials array
trialsarray = strcat('CTFeog = get(gcbf,''Userdata'');',...
  'CTFeog.TRIALS = str2num(get(CTFeog.handles.EctfTRIALS,''String'')); ',...
  'set(CTFeog.handles.EctfTRIALS,''Value'',CTFeog.TRIALS); ',...
  'set(gcbf,''Userdata'',CTFeog); clear CTFeog;');
G.EctfTRIALS = uicontrol('Parent',GUI,'Style','edit','Units','Normalized',Font, ...
  'Position',[.40 height .50 boxdepth], ...
  'String',num2str(CTFeog.TRIALS),'Value',CTFeog.TRIALS,...
  'TooltipString','Trials array', ...
  'Callback', trialsarray );
G.EctfTRIALSCLEAR = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized',Font, ...
  'Position',[.90 height .05 boxdepth], ...
  'String','Clear',...
  'TooltipString','Clear trials array', ...
  'Callback', strcat('CTFeog = get(gcbf,''Userdata'');',...
  'CTFeog.TRIALS = []; ',...
  'set(CTFeog.handles.EctfTRIALS,''String'',''''); ',...
  'set(CTFeog.handles.EctfTRIALS,''Value'',CTFeog.TRIALS); ',...
  'set(gcbf,''Userdata'',CTFeog); clear CTFeog;') );


% % PLOT: Load & plot the data!
% G.Bplot = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
%   'Position',[.20 .01 .18 .2],...
%   'String','PLOT','BusyAction','queue',...
%   'TooltipString','Plot the EEG data and return p struct.',...
%   'BackgroundColor',[0.0 0.5 0.0],...
%   'ForegroundColor', [1 1 1], 'HorizontalAlignment', 'center',...
%   'Callback',strcat('CTFeog = get(gcbf,''Userdata'');',...
%   'p = ctf_read_gui(CTFeog.p,''plot'');',...
%   'clear CTFeog;'));
% 
% % Save As
% G.Bsave = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
%   'Position',[.40 .01 .18 .2],'HorizontalAlignment', 'center',...
%   'String','SAVE AS','TooltipString','EEG File Conversion Tool (not implemented yet)',...
%   'BusyAction','queue',...
%   'Visible','off',...
%   'BackgroundColor',[0.0 0.0 0.75],...
%   'ForegroundColor', [1 1 1], 'HorizontalAlignment', 'center',...
%   'Callback',strcat('CTFeog = get(gcbf,''Userdata'');',...
%   'p = ctf_read_gui(CTFeog.p,''save'');',...
%   'clear CTFeog;'));

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
CTFeog.gui = GUI;          
CTFeog.handles = G;
CTFeog.ctf = ctf;
set(GUI,'Userdata',CTFeog);
set(GUI,'HandleVisibility','callback');

return








% ---------------------------------------------------------

function classTrials = extract_classTrials(CTFeog),
  
  ctf = CTFeog.ctf;
  
  ctf.class = ctf_read_classfile(ctf);
  
  if ctf.class.number,
    
    classNames = {ctf.class.data(:).name};
    classTrials = cell(1,ctf.class.number);
    for c = 1:ctf.class.number,
      classTrials{c} = ctf.class.data(c).trials;
    end
    
    % ------------------------------

    cd(ctf.folder)

    % -----------------------------
    % This section here is based on prior extraction of the STIM, EEG121 and
    % EEG122 channels, it could be replaced with the ctf_read_meg4 function
    % to do this directly

    stimeegFile = fullfile(ctf.folder,'WA_HiN_alltrials_1250msec_stimeeg.mat');
    stimeegData = load(stimeegFile);
    ctf.data = stimeegData.data;
    ctf.sensor.label = stimeegData.labels
    % -----------------------------

    stim = strmatch('STIM',ctf.sensor.label);
    veog = strmatch('EEG121',ctf.sensor.label);
    heog = strmatch('EEG122',ctf.sensor.label);

    blPoints = ctf.setup.pretrigger_samples;
    timeSec = ctf.setup.time_sec;
    epochs = size(ctf.data,3);


  end
  
  % For continuous data, use arbitrary epoch lengths for viewing
  %i = 1;
  %epochs(i,:) = 1:10000;
  %timeSec = epochs(i,:) * ctf.setup.sample_sec;
  %while (i * 10000) < size(ctf.data,1),
  %    epochs(i,:) = (i * 10000) + epochs(1,:);
  %    i = i + 1;
  %end


  return
