% pop_biosig() - import data files into EEGLAB using BIOSIG toolbox
%
% Usage:
%   >> OUTEEG = pop_biosig; % pop up window
%   >> OUTEEG = pop_biosig( filename, channels, type);
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
%  'blockepoch'  - ['on'|'off'] force importing continuous data. Default is 'on'
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
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, command, dat] = pop_biosig(filename, varargin); 
EEG = [];
command = '';

if ~plugin_askinstall('Biosig', 'sopen'), return; end;
biosigpathfirst;
    
if nargin < 1
	% ask user
	[filename, filepath] = uigetfile('*.*', 'Choose a data file -- pop_biosig()'); %%% this is incorrect in original version!!!!!!!!!!!!!!
    drawnow;
    
	if filename == 0 return; end;
	filename = [filepath filename];
    
    % look if MEG
    % -----------
    if length(filepath)>4
        if strcmpi(filepath(end-3:end-1), '.ds'), filename = filepath(1:end-1); end;
    end;
    
    % open file to get infos
    % ----------------------
    disp('Reading data file header...');
    dat = sopen(filename, 'r', [], 'OVERFLOWDETECTION:OFF');
    if ~isfield(dat, 'NRec')
        error('Unsuported data format');
    end;
    
    % special BIOSEMI
    % ---------------
    eeglab_options;
    if strcmpi(dat.TYPE, 'BDF')
        disp(upper('We highly recommend that you choose a reference channel IF these are Biosemi data'));
        disp(upper('(e.g., a mastoid or other channel). Otherwise the data will lose 40 dB of SNR!'));
        disp('For more information, see <a href="http://www.biosemi.com/faq/cms&drl.htm">http://www.biosemi.com/faq/cms&drl.htm</a>');
    end;
    uilist = { { 'style' 'text' 'String' 'Channel list (defaut all):' } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'String' [ 'Data range (in seconds) to read (default all [0 ' int2str(dat.NRec) '])' ] } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'String' 'Extract event' } ...
                 { 'style' 'checkbox' 'string' '' 'value' 1 'enable' 'on' } {} ...
                 { 'style' 'text' 'String' 'Import anotations (EDF+ only)' } ...
                 { 'style' 'checkbox' 'string' '' 'value' 1 'enable' 'on' } {} ...
                 { 'style' 'text' 'String' 'Force importing continuous data' 'value' 1} ...
                 { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
                 { 'style' 'text' 'String' 'Reference chan(s) indices - required for BIOSEMI' } ...
                 { 'style' 'edit' 'string' '' } ...
                  { 'style' 'checkbox' 'String' 'Import as memory mapped file (use if out of memory error)' 'value' option_memmapdata } };
    geom = { [3 1] [3 1] [3 0.35 0.5] [3 0.35 0.5] [3 0.35 0.5] [3 1] [1] };

    result = inputgui( geom, uilist, 'pophelp(''pop_biosig'')', ...
                                 'Load data using BIOSIG -- pop_biosig()');
    if length(result) == 0 return; end;
    
    % decode GUI params
    % -----------------
    options = {};
    if ~isempty(result{1}), options = { options{:} 'channels'   eval( [ '[' result{1} ']' ] ) }; end;
    if ~isempty(result{2}), options = { options{:} 'blockrange' eval( [ '[' result{2} ']' ] ) }; end;
    if length(result) > 2
        if ~result{3},          options = { options{:} 'importevent'  'off'  }; end;
        if ~result{4},          options = { options{:} 'importannot'  'off'  }; end;
        if  result{5},          options = { options{:} 'blockepoch'   'off' }; end;
        if ~isempty(result{6}), options = { options{:} 'ref'        eval( [ '[' result{6} ']' ] ) }; end;
        if  result{7},          options = { options{:} 'memorymapped' 'on' }; end;
    end;
    if length(eval( [ '[' result{6} ']' ] )) > 1
        options = { options{:} 'refoptions' { 'keepref' 'off' } };
    end;
else
    options = varargin;
end;

% decode imput parameters
% -----------------------
g = finputcheck( options, { 'blockrange'   'integer' [0 Inf]    [];
                            'channels'     'integer' [0 Inf]    [];
                            'ref'          'integer' [0 Inf]    [];
                            'refoptions'   'cell'    {}             { 'keepref' 'on' };
                            'rmeventchan'  'string'  { 'on';'off' } 'on';
                            'importevent'  'string'  { 'on';'off' } 'on';
                            'importannot'  'string'  { 'on';'off' } 'on';
                            'memorymapped' 'string'  { 'on';'off' } 'off';
                            'blockepoch'   'string'  { 'on';'off' } 'off' }, 'pop_biosig');
if isstr(g), error(g); end;

% import data
% -----------
EEG = eeg_emptyset;
[dat DAT interval] = readfile(filename, g.channels, g.blockrange, g.memorymapped);

if strcmpi(g.blockepoch, 'off')
    dat.NRec = 1;
end;
    
EEG = biosig2eeglab(dat, DAT, interval, g.channels, strcmpi(g.importevent, 'on'));

if strcmpi(g.rmeventchan, 'on') & strcmpi(dat.TYPE, 'BDF') & isfield(dat, 'BDF')
    if size(EEG.data,1) >= dat.BDF.Status.Channel, 
        disp('Removing event channel...');
        EEG.data(dat.BDF.Status.Channel,:) = []; 
        if ~isempty(EEG.chanlocs) && length(EEG.chanlocs) >= dat.BDF.Status.Channel
            EEG.chanlocs(dat.BDF.Status.Channel) = [];
        end;
    end;
    EEG.nbchan = size(EEG.data,1);
end;

% rerefencing
% -----------
if ~isempty(g.ref)
    disp('Re-referencing...');
    EEG = pop_reref(EEG, g.ref, g.refoptions{:});
%     EEG.data = EEG.data - repmat(mean(EEG.data(g.ref,:),1), [size(EEG.data,1) 1]);
%     if length(g.ref) == size(EEG.data,1)
%         EEG.ref  = 'averef';
%     end;
%     if length(g.ref) == 1
%         disp([ 'Warning: channel ' int2str(g.ref) ' is now zeroed (but still present in the data)' ]);
%     else
%         disp([ 'Warning: data matrix rank has decreased through re-referencing' ]);
%     end;
end;

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
            end;
        end;
    end;
end;

% convert data to single if necessary
% -----------------------------------
EEG = eeg_checkset(EEG,'makeur');   % Make EEG.urevent field
EEG = eeg_checkset(EEG);

% history
% -------
if isempty(options)
    command = sprintf('EEG = pop_biosig(''%s'');', filename); 
else
    command = sprintf('EEG = pop_biosig(''%s'', %s);', filename, vararg2str(options)); 
end;    

% Checking if str2double is on top of the path
biosigpathlast;

% ---------
% read data
% ---------
function [dat DAT interval] = readfile(filename, channels, blockrange, memmapdata);

if isempty(channels), channels = 0; end;
dat = sopen(filename, 'r', channels,'OVERFLOWDETECTION:OFF');

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
    
    if isempty(blockrange), blockrange = [0 dat.NRec]; end;
    blockrange(2) = min(blockrange(2), dat.NRec);
    allblocks = [blockrange(1):inc:blockrange(end)];
    count = 1;
    for bind = 1:length(allblocks)-1
        TMPDAT=sread(dat, (allblocks(bind+1)-allblocks(bind))*dat.Dur, allblocks(bind)*dat.Dur);
        if bind == 1
            DAT = mmo([], [size(TMPDAT,2) (allblocks(end)-allblocks(1))*dat.SPR]);
        end;
        DAT(:,count:count+length(TMPDAT)-1) = TMPDAT';
        count = count+length(TMPDAT);
    end;
    sclose(dat);
end;

if ~isempty(blockrange)
     interval(1) = blockrange(1) * dat.SampleRate(1) + 1;
     interval(2) = blockrange(2) * dat.SampleRate(1);
else interval = [];
end
