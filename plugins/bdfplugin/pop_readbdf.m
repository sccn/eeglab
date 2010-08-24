% pop_readbdf() - Read Biosemi 24-bit BDF file 
%
% Usage:
%   >> EEG = pop_readbdf;             % an interactive window pops up
%   >> EEG = pop_readbdf( filename ); % no pop-up window 
%   >> EEG = pop_readbdf( filename, range, eventchans, ref, delchan);
%
% Graphical interface:
%   "Data block range to read" - {edit box] see command line 
%                    parameter 'range' below.
%   "Event channel index(s)" - {edit box] see command line 
%                    parameter 'eventchans' below.
%   "Index(s) of reference channel(s)" - {edit box] see command line 
%                    parameter 'ref' below.
% Inputs:
%   filename       - Biosemi 24-bit BDF file name
%
% Optional input:
%   range          - [min max] integer range of data blocks to import.
%                    Default is empty -> import all data blocks.
%   eventchans     - [integer] event channel index(s). Default: none.
%   ref            - [integer] channel index or index(s) for the reference.
%                    Reference channels are not removed from the data,
%                    allowing easy re-referencing. If more than one
%                    channel, data are referenced to the average of the
%                    indexed channels. WARNING! Biosemi Active II data 
%                    are recorded reference-free, but LOSE 40 dB of SNR 
%                    if no reference is used!. If you do not know which
%                    channel to use, pick one and then re-reference after 
%                    the channel locations are read in. {default: none}
%   delchan        - ['on'|'off'] delete event channel. Default if 'on'\
%
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 13 March 2002
%
% See also: openbdf(), readbdf()

% Copyright (C) 13 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%adbdf( filename ); % no pop-up window 
%   >> EEG = pop_readbdf( filename, range, eventchans, ref, delchan);
%
% Graphical interface:
%   "Data block range to read" - {edit box] see command line 
%                    parameter 'range' below.
%   "Event channel index(s)" - {edit box] see command line 
%                    parameter 'eventchans' below.
%   "Index(s) of reference channel(s)" - {edit box] see command line 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% programmed from pop_readedf() version 1.15

function [EEG, command] = pop_readbdf(filename, blockrange, eventchans, ref, delchan); 
EEG = [];
command = '';

if nargin < 2
    blockrange = [];
end;
if nargin < 3
    eventchans = [];
end;
if nargin < 4
    ref = [];
end;
if nargin < 5
    delchan = 'on';
end;

if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.BDF;*.bdf', 'Choose an BDF file -- pop_readbdf()'); 
    drawnow;
    
	if filename == 0 return; end;
	filename = [filepath filename];
    disp('We highly recommend that you choose a reference channel IF these are Biosemi data');
    disp('(e.g., a mastoid or other channel). Otherwise the data will lose 40 dB of SNR!');
    
    dat = openbdf(filename);
    promptstr    = { [ 'Data block range to read (default all [1 ' int2str(dat.Head.NRec) '])' ]
                     [ 'Event channel number(s) (default:none [1-' int2str(dat.Head.NS)  '])' ]
                     'If Biosemi data, reference chan(s) number(s)' 
                     'Delete event channel' };
    inistr       = { '' int2str(dat.Head.NS) '' 'on' };
    result       = inputdlg2( promptstr, 'Import BDF file -- pop_readbdf()', 1,  inistr, 'pop_readbdf');
    if length(result) == 0 return; end;
    blockrange   = eval( [ '[' result{1} ']' ] );
    eventchans   = eval( [ '[' result{2} ']' ] );
    ref          = eval( [ '[' result{3} ']' ] );
    delchan      = result{4};
end;

% load datas
% ----------
EEG = eeg_emptyset;
fprintf('Reading BDF data in 24-bit format...\n');
dat = openbdf(filename);
%dat.Head = sopen(filename);
if isempty(dat.Head.NRec)
    dat.Head.NRec = 100;
end;
if isempty(blockrange)
    blockrange = [1 dat.Head.NRec];
end;
vectrange = [blockrange(1):min(blockrange(2), dat.Head.NRec)];
DAT=readbdf(dat,vectrange);
EEG.nbchan          = size(DAT.Record,1);
EEG.srate           = dat.Head.SampleRate(1);
EEG.data            = DAT.Record;
EEG.pnts            = size(EEG.data,2);
EEG.trials          = 1;
EEG.setname 		= 'BDF file';
disp('Event information might be encoded in the last channel');
disp('To extract these events, use menu File > Import event info > From data channel'); 
EEG.filename        = filename;
EEG.filepath        = '';
EEG.xmin            = 0; 
EEG.chanlocs        = struct('labels', cellstr(dat.Head.Label));
EEG = eeg_checkset(EEG);

% extract events
% --------------
disp('Extracting events...');
warning off;
if ~isempty(eventchans)
    if any(eventchans > EEG.nbchan) | any(eventchans < 1)
        error('Event channel index or indices out of range');
    end;
    disp('Applying 8-bit masks to event channels (only for BIOSEMI files).');
    %EEG.data(eventchans,:) = bitand(uint32(EEG.data(eventchans,:)), 32767);
    
    %thiscode = 0;lastcode=0;
    %for p = 1:size(EEG.data,2)-1
    %    prevcode = thiscode;
    %    thiscode = mod(EEG.data(end,p),256*256);   % andrey's code - 16 bits 
    %    if (thiscode ~= 0) && (thiscode~=prevcode) && (thiscode~=lastcode) % fix to avoid repeated codes (per Ying's demand)
    %        EEG.event(end+1).latency =  p;
    %        EEG.event(end).type = thiscode;
    %        lastcode = thiscode;
    %    end;
    %end;
    thiscode = 0;
    for p = 1:size(EEG.data,2)-1
        prevcode = thiscode;
        thiscode = mod(EEG.data(end,p),256*256);   % andrey's code - 16 bits 
        if (thiscode ~= 0) && (thiscode~=prevcode) 
            EEG.event(end+1).latency =  p;
            EEG.event(end).type = thiscode;
        end;
    end;
    if strcmpi(delchan, 'on')
        EEG = pop_select(EEG, 'nochannel', size(EEG.data,1));
    end;
    EEG = eeg_checkset(EEG, 'makeur');
    
    %EEG = pop_chanevent(EEG, eventchans, 'edge', 'leading', 'delchan', delchan);
end;
warning on;

% rerefencing
% -----------
if ~isempty(ref)
    disp('Re-referencing...');
    if ~isempty(eventchans) & any(ref > eventchans(1))
        error('Event channel index cannot be before ref channel index, contact arno@salk.edu');
    end;
    EEG.data = EEG.data - repmat(mean(EEG.data(ref,:),1), [size(EEG.data,1) 1]);
    if length(ref) == 1
        disp([ 'Warning: channel ' int2str(ref) ' is now zeroed (but still present in the data)' ]);
    else
        disp([ 'Warning: data matrix rank has decreased through re-referencing' ]);
    end;
end;

% history
% -------
if ~isempty(blockrange) | ~isempty(eventchans)
    command = sprintf('EEG = pop_readbdf(''%s'', %s);', filename, vararg2str({blockrange eventchans ref})); 
else
    command = sprintf('EEG = pop_readbdf(''%s'');', filename); 
end;    

return;
