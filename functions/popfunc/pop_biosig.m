% pop_biosig() - import data files into EEGLAB using BIOSIG toolbox
%
% Usage:
%   >> OUTEEG = pop_biosig; % pop up window
%   >> OUTEEG = pop_biosig( filename, 'key', val, ...);
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'channels'   - [integer array] list of channel indices
%   'blockrange' - [min max] integer range of data blocks to import, in seconds.
%                  Entering [0 3] will import the first three blocks of data.
%                  Default is empty -> import all data blocks. 
%  'importevent' - ['on'|'off'] import events. Default if 'on'.
%  'importannot' - ['on'|'off'] import annotations (EDF+ only). Default if 'on'
%  'importmex'   - ['on'|'off'] import events with Biosig mexSLOAD as an alternative. Default if 'off'
%  'overflow'    - ['on'|'off'] overflow detection. Default is 'off'
%  'uncalibrated' - ['on'|'off'] import uncalibrated data. Default is 'off'
%  'blockepoch'  - ['on'|'off'] force importing continuous data. Default is 'on'
%  'bdfeventmode' - [interger] see bdf2biosig_events function help. Default is 4.
%  'ref'         - [integer] channel index or index(s) for the reference.
%                  Reference channels are not removed from the data,
%                  allowing easy re-referencing. If more than one
%                  channel, data are referenced to the average of the
%                  indexed channels. WARNING! Biosemi Active II data 
%                  are recorded reference-free, but LOSE 40 dB of SNR 
%                  if no reference is used!. If you do not know which
%                  channel to use, pick one and then re-reference after 
%                  the channel locations are read in. {default: none}.
%                  For more information see http://www.biosemi.com/faq/cms&drl.htm
%  'refoptions'  - [Cell] Option for the pop_reref function. Default is to 
%                  remove the reference channel if there is one of them and to 
%                  keep it if there are several of them from the graphic
%                  interface. From the command line default option is to 
%                  keep the reference channel.
%  'rmeventchan' - ['on'|'off'] remove event channel after event 
%                  extraction. Default is 'on'.
%  'memorymapped' - ['on'|'off'] import memory mapped file (useful if 
%                  encountering memory errors). Default is 'off'.
%
% Outputs:
%   OUTEEG   - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Oct. 29, 2003-
%
% Note: BIOSIG toolbox must be installed. Download BIOSIG at 
%       http://biosig.sourceforge.net
%       Contact a.schloegl@ieee.org for troubleshooting using BIOSIG.

% Copyright (C) 2003 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [ALLEEG, command, dat] = pop_biosig(filename, varargin)
ALLEEG = [];
command = '';

if ~plugin_askinstall('Biosig', 'sopen'), return; end
biosigpathfirst;
    
saveData = false;
if nargin < 1
	% ask user
	[filename, filepath] = uigetfile('*.*', 'Choose a data file -- pop_biosig()', 'multiselect', 'on'); %%% this is incorrect in original version!!!!!!!!!!!!!!
    drawnow;
    
	if isequal(filename,0) return; end
    
    if iscell(filename)
        buttonName = questdlg2([ 'Do you want to automatically save imported datasets?' 10 ...
            '(the name will remain the same as the original dataset' 10 ...
            'and the .set extension will be used)' ], 'pop_biosig() - import dataset(s)', 'Cancel', 'No thanks', 'Yes Save', 'Yes Save');
        switch buttonName
            case 'Cancel', return;
            case 'No thanks', saveData = false;
            otherwise saveData = true;
        end
    else
        filename = { filename };
    end

    % look if MEG
    % -----------
    if length(filepath)>4
        if strcmpi(filepath(end-3:end-1), '.ds')
            filename = { filepath(1:end-1) }; % should not be able to select more than one MEG file anyway
        end
    end
    for iFile = 1:length(filename)
        filename{iFile} = fullfile(filepath, filename{iFile});
    end
    
    % open file to get infos
    % ----------------------
    disp('Reading data file header...');
    dat = sopen(filename{1}, 'r', [], 'OVERFLOWDETECTION:OFF');
    if ~isfield(dat, 'NRec')
        error('Unsuported data format');
    end
    
    % special BIOSEMI
    % ---------------
    eeglab_options;
    if strcmpi(dat.TYPE, 'BDF')
        disp(upper('We highly recommend that you choose a reference channel IF these are Biosemi data'));
        disp(upper('(e.g., a mastoid or other channel). Otherwise the data will lose 40 dB of SNR!'));
        disp('For more information, see <a href="http://www.biosemi.com/faq/cms&drl.htm">http://www.biosemi.com/faq/cms&drl.htm</a>');
    end
    checkmex = [ 'if ~exist(''mexSLOAD''), set(gcbo, ''value'', 0); ' ...
                 'warndlg2([ ''mexSLOAD not found in path. It needs to be installed.'' 10 ' ...
                 '''It is easier to use this option on Windows where mexSLOAD'' 10 ' ...
                 '''is automatically installed with BIOSIG.'' ]); end' ];
    uilist = { { 'style' 'text' 'String' 'Channel list (default all):' } ...
                 { 'style' 'edit' 'string' '' 'tag' 'channels' } ...
                 { 'style' 'text' 'String' [ 'Data range (in seconds) to read (default all [0 ' int2str(dat.NRec) '])' ] } ...
                 { 'style' 'edit' 'string' '' 'tag' 'blockrange'  } ...
                 { 'style' 'text' 'String' 'Reference chan(s) indices - required for BIOSEMI' } ...
                 { 'style' 'edit' 'string' ''  'tag' 'ref' } ...
                 { 'style' 'text' 'String' 'Event mode import (see bdf2biosig_events function help)' } ...
                 { 'style' 'edit' 'string' '4'  'tag' 'bdfmode' } ...         
                 { 'style' 'checkbox' 'string' 'Extract event from marker channel' 'value' 1 'tag' 'importevent' } ...
                 { 'style' 'checkbox' 'String' 'Force continuous data when data is stored in blocks' 'value' 1 'tag' 'blockepoch' } ...
                 { 'style' 'checkbox' 'String' 'Import as memory mapped file (use if out of memory error)' 'value' option_memmapdata 'tag' 'memorymapped' } ...
                 { 'style' 'checkbox' 'String' 'Enable overflow detection' 'value' 0 'tag' 'overflow' } ...
                 { 'style' 'checkbox' 'String' 'Import uncalibrated data' 'value' 0 'tag' 'uncalibrated' } ...
                 { 'style' 'checkbox' 'string' 'Import EDF+ anotations (try also mexSLOAD below)' 'value' 1 'enable' 'on' 'tag' 'importannot' } ...
                 { 'style' 'checkbox' 'string' 'Import using alternative BIOSIG method (mexSLOAD)' 'value' 0 'callback' checkmex 'tag' 'importmex' } };
    geom = { [3 1] [3 1] [3 1] [3 1] [1] [1] [1] [1] [1] [1] [1] };

    [~,~,~,result] = inputgui( geom, uilist, 'pophelp(''pop_biosig'')', ...
                                 'Load data using BIOSIG -- pop_biosig()');
    if length(result) == 0 return; end
    
    % decode GUI params
    % -----------------
    options = {};
    if ~isempty(result.channels)  ,   options = { options{:} 'channels'     eval( [ '[' result.channels ']' ] ) }; end
    if ~isempty(result.blockrange),   options = { options{:} 'blockrange'   eval( [ '[' result.blockrange ']' ] ) }; end
    if ~isempty(result.ref       ),   options = { options{:} 'ref'          eval( [ '[' result.ref ']' ] ) }; end
    if ~isequal(result.bdfmode, '4'), options = { options{:} 'bdfeventmode' eval( [ '[' result.bdfmode ']' ] ) }; end
    
    % default to 1
    if ~result.importevent,    options = { options{:} 'importevent'  'off'  }; end
    if ~result.blockepoch ,    options = { options{:} 'blockepoch'   'off'  }; end
    if ~result.importannot ,   options = { options{:} 'importannot'  'off'  }; end
    
    % default to 0
    if result.memorymapped,    options = { options{:} 'memorymapped'  'on'  }; end
    if result.importmex   ,    options = { options{:} 'importmex'     'on'  }; end
    if result.overflow    ,    options = { options{:} 'overflow'      'on'  }; end
    if result.uncalibrated,    options = { options{:} 'uncalibrated'  'on'  }; end
    
    if length(eval( [ '[' result.ref ']' ] )) > 1
        options = { options{:} 'refoptions' { 'keepref' 'off' } };
    end
else
    options = varargin;
end

% decode input parameters
% -----------------------
g = finputcheck( options, { 'blockrange'   'integer' [0 Inf]    [];
                            'channels'     'integer' [0 Inf]    [];
                            'ref'          'integer' [0 Inf]    [];
                            'bdfeventmode' 'integer' [0 Inf]    4;
                            'refoptions'   'cell'    {}             { 'keepref' 'on' };
                            'rmeventchan'  'string'  { 'on';'off' } 'on';
                            'overflow'     'string'  { 'on';'off' } 'off';
                            'uncalibrated' 'string'  { 'on';'off' } 'off';
                            'importevent'  'string'  { 'on';'off' } 'on';
                            'importannot'  'string'  { 'on';'off' } 'on';
                            'importmex'   'string'  { 'on';'off' }  'off';
                            'memorymapped' 'string'  { 'on';'off' } 'off';
                            'blockepoch'   'string'  { 'on';'off' } 'off' }, 'pop_biosig');
if ischar(g), error(g); end

if ~iscell(filename) filename = { filename }; end

for iFile = 1:length(filename)
    % import data
    % -----------
    EEG = eeg_emptyset;
    [dat, DAT, interval] = readfile(filename{iFile}, g.channels, g.blockrange, g.memorymapped, g.bdfeventmode, g.overflow, g.uncalibrated);
    
    if strcmpi(g.blockepoch, 'off')
        dat.NRec = 1;
    end
    
    EEG = biosig2eeglab(dat, DAT, interval, g.channels, strcmpi(g.importevent, 'on'), strcmpi(g.importannot, 'on'));
    
    if strcmpi(g.rmeventchan, 'on') && strcmpi(dat.TYPE, 'BDF') && isfield(dat, 'BDF')
        if size(EEG.data,1) >= dat.BDF.Status.Channel
            disp('Removing event channel...');
            EEG.data(dat.BDF.Status.Channel,:) = [];
            if ~isempty(EEG.chanlocs) && length(EEG.chanlocs) >= dat.BDF.Status.Channel
                EEG.chanlocs(dat.BDF.Status.Channel) = [];
            end
        end
        EEG.nbchan = size(EEG.data,1);
    end
    
    % rerefencing
    % -----------
    if ~isempty(g.ref)
        disp('Re-referencing...');
        EEG = pop_reref(EEG, g.ref, g.refoptions{:});
    end
    
    % test if annotation channel is present
    % -------------------------------------
    if isfield(dat, 'EDFplus') && strcmpi(g.importannot, 'on')
        tmpfields = fieldnames(dat.EDFplus);
        for ind = 1:length(tmpfields)
            tmpdat = getfield(dat.EDFplus, tmpfields{ind});
            if length(tmpdat) == EEG.pnts
                EEG.data(end+1,:) = tmpdat;
                EEG.nbchan        = EEG.nbchan+1;
                if ~isempty(EEG.chanlocs)
                    EEG.chanlocs(end+1).labels = tmpfields{ind};
                end
            end
        end
    end
    
    % import using Biosig mexSLOAD method (Cedric edits 2/23/2021)
    if strcmpi(g.importmex, 'on')
        [s,HDR] = mexSLOAD(filename);
        
        %Get correct event names contained in CodeDesc
        num_ev_type = unique(HDR.EVENT.TYP);
        num_ev_name = unique(HDR.EVENT.CodeDesc);
        if ~isempty(HDR.EVENT.CodeDesc) && length(num_ev_type) == length(num_ev_name)
            for iEvent = 1:length(HDR.EVENT.TYP)
                EEG.event(iEvent).type = char(HDR.EVENT.CodeDesc(HDR.EVENT.TYP(iEvent)));
                EEG.event(iEvent).latency = HDR.EVENT.POS(iEvent);
            end
        else
            for iEvent = 1:length(HDR.EVENT.TYP)
                EEG.event(iEvent).type = HDR.EVENT.TYP(iEvent);
                EEG.event(iEvent).latency = HDR.EVENT.POS(iEvent);
                EEG.event(iEvent).urevent = iEvent;
            end
            warning('Inconsistency between event types and names or event names were not found.');
            warning('Check your events or convert your data format with EDFBrowser: https://www.teuniz.net/edfbrowser/');
        end
    end
    
    % check and store data
    % --------------------
%     EEG = eeg_checkset(EEG, 'makeur');   % Make EEG.urevent field
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG = eeg_checkset(EEG);
    if saveData
        EEG = pop_saveset(EEG, [ filename{iFile}(1:end-4) '.set' ]);
        EEG(1).saved = 'justloaded';
    end
    if iFile == 1
        ALLEEG = EEG;
    else
        ALLEEG(iFile) = EEG;
    end
    
end

% history
% -------
if length(filename) == 1
    command = sprintf('EEG = pop_biosig(''%s''', filename{1});
else
    command = sprintf('EEG = pop_biosig(%s', vararg2str(filename));
end
if isempty(options)
    command = [ command ');' ];
else
    command = [ command sprintf(', %s);', vararg2str(options)) ];
end

% Checking if str2double is on top of the path
biosigpathlast;

% ---------
% read data
% ---------
function [dat, DAT, interval] = readfile(filename, channels, blockrange, memmapdata, bdfeventmode, overflow, uncalibrated)

if isempty(channels), channels = 0; end
strmode = '';
if strcmpi(overflow, 'off')
    strmode = 'OVERFLOWDETECTION:OFF';
end
if strcmpi(uncalibrated, 'on')
    if ~isempty(strmode) strmode = [ strmode ';' ]; end
    strmode = [ strmode 'UCAL'];
end
if ~isequal(bdfeventmode, 4)
    if ~isempty(strmode) strmode = [ strmode ';' ]; end
    strmode = [ strmode 'BDF:[' num2str(bdfeventmode) ']' ];
end
fprintf('sopen mode is "%s"\n', strmode);
dat = sopen(filename, 'r', channels, strmode);

if strcmpi(memmapdata, 'off')
    fprintf('Reading data in %s format...\n', dat.TYPE);

    if ~isempty(blockrange)
        newblockrange    = blockrange;
%         newblockrange    = newblockrange*dat.Dur;    
        DAT=sread(dat, newblockrange(2)-newblockrange(1), newblockrange(1));
    else 
        DAT=sread(dat, Inf);% this isn't transposed in original!!!!!!!!
        newblockrange    = [];
    end
    sclose(dat);
else
    fprintf('Reading data in %s format (file will be mapped to memory so this may take a while)...\n', dat.TYPE);
    inc = ceil(250000/(dat.NS*dat.SPR)); % 1Mb block
    
    if isempty(blockrange), blockrange = [0 dat.NRec]; end
    blockrange(2) = min(blockrange(2), dat.NRec);
    allblocks = [blockrange(1):inc:blockrange(end)];
    count = 1;
    for bind = 1:length(allblocks)-1
        TMPDAT=sread(dat, (allblocks(bind+1)-allblocks(bind))*dat.Dur, allblocks(bind)*dat.Dur);
        if bind == 1
            DAT = mmo([], [size(TMPDAT,2) (allblocks(end)-allblocks(1))*dat.SPR]);
        end
        DAT(:,count:count+length(TMPDAT)-1) = TMPDAT';
        count = count+length(TMPDAT);
    end
    sclose(dat);
end

if ~isempty(blockrange)
     interval(1) = blockrange(1) * dat.SampleRate(1) + 1;
     interval(2) = blockrange(2) * dat.SampleRate(1);
else interval = [];
end
