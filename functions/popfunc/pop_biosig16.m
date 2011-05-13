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
%   'ref'        - [integer] channel index or index(s) for the reference.
%                  Reference channels are not removed from the data,
%                  allowing easy re-referencing. If more than one
%                  channel, data are referenced to the average of the
%                  indexed channels. WARNING! Biosemi Active II data 
%                  are recorded reference-free, but LOSE 40 dB of SNR 
%                  if no reference is used!. If you do not know which
%                  channel to use, pick one and then re-reference after 
%                  the channel locations are read in. {default: none}
%   'rmeventchan' - ['on'|'off'] remove event channel after event 
%                  extraction. Default is 'on'.
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

function [EEG, command] = my_pop_biosig(filename, varargin); 
EEG = [];
command = '';

if nargin < 1
	% ask user
	[filename, filepath] = uigetfile('*.*', 'Choose an BDF file -- pop_biosig()'); %%% this is incorrect in original version!!!!!!!!!!!!!!
    drawnow;
    
	if filename == 0 return; end;
	filename = [filepath filename];
    
    % open file to get infos
    % ----------------------
    disp('Reading data file header...');
    dat = sopen(filename);
    
    % special BIOSEMI
    % ---------------
    if strcmpi(dat.TYPE, 'BDF')
        disp('We highly recommend that you choose a reference channel IF these are Biosemi data');
        disp('(e.g., a mastoid or other channel). Otherwise the data will lose 40 dB of SNR!');
    end;
    uilist = { { 'style' 'text' 'String' 'Channel list (defaut all):' } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'String' [ 'Data range (in seconds) to read (default all [0 ' int2str(dat.NRec) '])' ] } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'String' 'Extract event - cannot be unset (set=yes)' } ...
                 { 'style' 'checkbox' 'string' '' 'value' 1 'enable' 'on' } {} ...
                 { 'style' 'text' 'String' 'Import continuous data (set=yes)' 'value' 1} ...
                 { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
                 { 'style' 'text' 'String' 'Reference chan(s) indices - required for BIOSEMI' } ...
                 { 'style' 'edit' 'string' '' } };
    geom = { [3 1] [3 1] [3 0.35 0.5] [3 0.35 0.5] [3 1] };
    result = inputgui( geom, uilist, 'pophelp(''pop_biosig'')', ...
                                 'Load data using BIOSIG -- pop_biosig()');
    if length(result) == 0 return; end;
    
    % decode GUI params
    % -----------------
    options = {};
    if ~isempty(result{1}), options = { options{:} 'channels'   eval( [ '[' result{1} ']' ] ) }; end;
    if ~isempty(result{2}), options = { options{:} 'blockrange' eval( [ '[' result{2} ']' ] ) }; end;
    if length(result) > 2
        if ~isempty(result{4}), options = { options{:} 'ref'        eval( [ '[' result{4} ']' ] ) }; end;
        if ~result{3},          options = { options{:} 'rmeventchan' 'off' }; end;
    end;
else
    options = varargin;
end;

% decode imput parameters
% -----------------------
g = finputcheck( options, { 'blockrange'  'integer' [0 Inf]    [];
                            'channels'    'integer' [0 Inf]    [];
                            'ref'         'integer' [0 Inf]    [];
                            'rmeventchan' 'string'  { 'on';'off' } 'on' }, 'pop_biosig');
if isstr(g), error(g); end;

% import data
% -----------
EEG = eeg_emptyset;
if ~isempty(g.channels)
     dat = sopen(filename, 'r', g.channels,'OVERFLOWDETECTION:OFF');
else dat = sopen(filename, 'r', 0,'OVERFLOWDETECTION:OFF');
end
fprintf('Reading data in %s format...\n', dat.TYPE);

if ~isempty(g.blockrange)
    newblockrange    = g.blockrange;
    newblockrange(2) = min(newblockrange(2), dat.NRec);
    newblockrange    = newblockrange*dat.Dur;    
    DAT=sread(dat, newblockrange(2)-newblockrange(1), newblockrange(1))';
else 
    DAT=sread(dat, Inf)';% this isn't transposed in original!!!!!!!!
    newblockrange    = [];
end
dat = sclose(dat);

% convert to seconds for sread
% ----------------------------
EEG.nbchan          = size(DAT,1);
EEG.srate           = dat.SampleRate(1);
EEG.data            = DAT; 
clear DAT;
% $$$ try  % why would you do the following???????  JO
% $$$     EEG.data            = EEG.data';
% $$$ catch,
% $$$     pack;
% $$$     EEG.data            = EEG.data';
% $$$ end;    
EEG.setname 		= sprintf('%s file', dat.TYPE);
EEG.comments        = [ 'Original file: ' filename ];
EEG.xmin            = 0; 
if strcmpi(dat.TYPE, 'BDF') || strcmpi(dat.TYPE, 'EDF')
    EEG.trials   = 1;
    EEG.pnts     = size(EEG.data,2);
else
    EEG.trials   = dat.NRec;
    EEG.pnts     = size(EEG.data,2)/dat.NRec;
end
if isfield(dat, 'Label') & ~isempty(dat.Label)
    EEG.chanlocs = struct('labels', cellstr(char(dat.Label)));
end
EEG = eeg_checkset(EEG);

% extract events % this part I totally revamped to work...  JO
% --------------
disp('Extracting events from last EEG channel...');
EEG.event = [];

% $$$ startval = mode(EEG.data(end,:)); % my code
% $$$ for p = 2:size(EEG.data,2)-1
% $$$     [codeout] = code(EEG.data(end,p));
% $$$     if EEG.data(end,p) > EEG.data(end,p-1) & EEG.data(end,p) >= EEG.data(end,p+1)
% $$$         EEG.event(end+1).latency =  p;
% $$$         EEG.event(end).type = bitand(double(EEG.data(end,p)-startval),255);
% $$$     end;
% $$$ end;

% lastout = mod(EEG.data(end,1),256);newevs = []; % andrey's code 8 bits
% codeout = mod(EEG.data(end,2),256);
% for p = 2:size(EEG.data,2)-1
%     nextcode = mod(EEG.data(end,p+1),256);
%     if codeout > lastout & codeout >= nextcode
%         newevs = [newevs codeout];
%         EEG.event(end+1).latency =  p;
%         EEG.event(end).type = codeout;
%     end;
%     lastout = codeout;
%     codeout = nextcode;
% end;

%lastout = mod(EEG.data(end,1),256*256);newevs = []; % andrey's code 16 bits
%codeout = mod(EEG.data(end,2),256*256);
%for p = 2:size(EEG.data,2)-1
%    nextcode = mod(EEG.data(end,p+1),256*256);
%    if (codeout > lastout) & (codeout >= nextcode)
%        newevs = [newevs codeout];
%        EEG.event(end+1).latency =  p;
%        EEG.event(end).type = codeout;
%    end;
%    lastout = codeout;
%    codeout = nextcode;
%end;

% Modifieded by Andrey (Aug.5,2008) to detect all non-zero codes: 
thiscode = 0;
for p = 1:size(EEG.data,2)-1
    prevcode = thiscode;
    thiscode = mod(EEG.data(end,p),256*256);   % andrey's code - 16 bits 
    if (thiscode ~= 0) && (thiscode~=prevcode) 
        EEG.event(end+1).latency =  p;
        EEG.event(end).type = thiscode;
    end;
end;

if strcmpi(g.rmeventchan, 'on')
    EEG.data(dat.BDF.Status.Channel,:) = [];
    EEG.nbchan = size(EEG.data,1);
    if ~isempty(EEG.chanlocs)
        EEG.chanlocs(dat.BDF.Status.Channel,:) = [];
    end;
end;
EEG = eeg_checkset(EEG, 'eventconsistency');

% $$$ if ~isempty(dat.EVENT)    
% $$$     if isfield(dat, 'out') % Alois fix for event interval does not work
% $$$         if isfield(dat.out, 'EVENT')
% $$$             dat.EVENT = dat.out.EVENT;
% $$$         end;
% $$$     end;
% $$$     if ~isempty(newblockrange)
% $$$         interval(1) = newblockrange(1) * dat.SampleRate(1) + 1;
% $$$         interval(2) = newblockrange(2) * dat.SampleRate(1);
% $$$     else interval = [];
% $$$     end
% $$$     EEG.event = biosig2eeglabevent(dat.EVENT, interval); % Toby's fix
% $$$     if strcmpi(g.rmeventchan, 'on') & strcmpi(dat.TYPE, 'BDF') & isfield(dat, 'BDF')
% $$$         disp('Removing event channel...');
% $$$         EEG.data(dat.BDF.Status.Channel,:) = [];
% $$$         EEG.nbchan = size(EEG.data,1);
% $$$         if ~isempty(EEG.chanlocs)
% $$$             EEG.chanlocs(dat.BDF.Status.Channel,:) = [];
% $$$         end;
% $$$     end;
% $$$     EEG = eeg_checkset(EEG, 'eventconsistency');
% $$$ else 
% $$$     disp('Warning: no event found. Events might be embeded in a data channel.');
% $$$     disp('         To extract events, use menu File > Import Event Info > From data channel');
% $$$ end;

% rerefencing
% -----------
if ~isempty(g.ref)
    disp('Re-referencing...');
    EEG.data = EEG.data - repmat(mean(EEG.data(g.ref,:),1), [size(EEG.data,1) 1]);
    if length(g.ref) == size(EEG.data,1)
        EEG.ref  = 'averef';
    end;
    if length(g.ref) == 1
        disp([ 'Warning: channel ' int2str(g.ref) ' is now zeroed (but still present in the data)' ]);
    else
        disp([ 'Warning: data matrix rank has decreased through re-referencing' ]);
    end;
end;

% convert data to single if necessary
% -----------------------------------
EEG = eeg_checkset(EEG,'makeur');   % Make EEG.urevent field
EEG = eeg_checkset(EEG);

% history
% -------
if isempty(options)
    command = sprintf('EEG = my_pop_biosig(''%s'');', filename); 
else
    command = sprintf('EEG = my_pop_biosig(''%s'', %s);', filename, vararg2str(options)); 
end;    
