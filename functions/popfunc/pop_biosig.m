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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.35  2009/05/07 20:29:33  arno
% edit message
%
% Revision 1.34  2009/05/07 00:01:23  arno
% now calling biosig2eeglab
%
% Revision 1.33  2009/05/01 01:06:32  arno
% Allow writing EDF/BDF/GDF files
%
% Revision 1.32  2009/04/21 21:24:45  arno
% Biosig fix for EDF
%
% Revision 1.31  2009/03/20 21:58:00  arno
% Put back old event format
%
% Revision 1.30  2009/02/12 21:55:24  arno
% fix new event channel etc...
%
% Revision 1.29  2008/11/13 00:41:49  arno
% close file later on
%
% Revision 1.28  2008/08/11 22:04:15  nima
% Andrey fixed a subtle bug that prevented events from being detected when they were in adjacent samples (not separated by zeros). Also, the function now can import 16 bit codes (1 to 65535) from the event channel.
%
% Revision 1.27  2008/06/30 22:19:03  arno
% julie's changes
%
% Revision 1.26  2008/05/08 22:47:27  elizabeth
% Test for empty chanlocs - Arno
%
% Revision 1.25  2008/04/17 18:33:17  nima
% replaced by arno files version of the file
%
% Revision 1.23  2007/08/14 20:07:29  arno
% fix for new BIOSIG version
%
% Revision 1.22  2007/08/02 22:45:39  arno
% back to previous
%
% Revision 1.21  2007/08/02 22:35:37  arno
% adding channel
% for statust field (events)
%
% Revision 1.20  2007/03/20 16:33:35  arno
% put back interval
%
% Revision 1.19  2007/03/20 16:28:32  arno
% remove toby's fix
%
% Revision 1.18  2007/03/20 02:05:03  arno
% checking dat.out.EVENT
%
% Revision 1.17  2007/03/08 03:08:02  toby
% Bug 216
%
% Revision 1.16  2006/12/16 07:36:21  toby
% Fixed partial data event problem, documentation update.
%
% Revision 1.15  2006/12/16 04:16:37  toby
% Help edit.
%
% Revision 1.14  2006/07/05 23:03:47  arno
% checking for field BDF
%
% Revision 1.13  2006/04/19 14:54:21  arno
% EDF continous files
%
% Revision 1.12  2006/01/25 21:11:29  arno
% fixing 2 typos3
%
% Revision 1.11  2006/01/25 00:38:08  arno
% fixing sclose
%
% Revision 1.10  2006/01/17 00:56:19  arno
% closing data file after reading
%
% Revision 1.9  2006/01/17 00:54:13  arno
% disable overflow detection3
%
% Revision 1.8  2006/01/17 00:47:08  arno
% 0verflow detection off
%
% Revision 1.7  2006/01/13 23:07:54  arno
% new generic function for all BIOSIG
%
% Revision 1.6  2006/01/13 22:23:01  arno
% same
%
% Revision 1.5  2006/01/13 22:21:32  arno
% special handling of BDF files
%
% Revision 1.4  2006/01/13 22:01:01  arno
% now uses biosig2eeglabevent
%
% Revision 1.3  2005/10/27 05:24:55  arno
% filename3
%
% Revision 1.2  2005/03/23 16:01:17  arno
% implmenting DUR and CHN
%
% Revision 1.1  2004/11/12 18:22:11  arno
% Initial revision
%
% Revision 1.1  2004/09/12 02:03:47  arnodelorme
% Adding EEGLAB folder with EEGLAB interface files
%
% Revision 1.4  2004/08/31 21:08:43  arno
% new messages
%
% Revision 1.3  2003/12/19 17:33:23  arno
% message to import data
%
% Revision 1.2  2003/12/19 17:28:50  arno
% importing events
%
% Revision 1.1  2003/12/19 17:18:43  arno
% Initial revision
%
% Revision 1.3  2003/10/29 18:53:31  arno
% text typos
%
% Revision 1.2  2003/10/29 18:49:31  arno
% debuging type
%
% Revision 1.1  2003/10/29 18:17:26  arno
% Initial revision
%

function [EEG, command, dat] = my_pop_biosig(filename, varargin); 
EEG = [];
command = '';

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
    dat = sopen(filename);
    uilist   = { { 'style' 'text' 'String' 'Channel list (defaut all):' } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'String' [ 'Data range (in seconds) to read (default all [0 ' int2str(dat.NRec) '])' ] } ...
                 { 'style' 'edit' 'string' '' }  };
    geom = { [3 1] [3 1] };
    
    % special BIOSEMI
    % ---------------
    if strcmpi(dat.TYPE, 'BDF')
        disp('We highly recommend that you choose a reference channel IF these are Biosemi data');
        disp('(e.g., a mastoid or other channel). Otherwise the data will lose 40 dB of SNR!');
    end;
    uilist = { uilist{:} ...
                 { 'style' 'text' 'String' 'Extract event - cannot be unset (set=yes)' } ...
                 { 'style' 'checkbox' 'string' '' 'value' 1 'enable' 'off' } {} ...
                 { 'style' 'text' 'String' 'Import continuous data (set=yes)' 'value' 1} ...
                 { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
                 { 'style' 'text' 'String' 'Reference chan(s) indices - required for BIOSEMI' } ...
                 { 'style' 'edit' 'string' '' } };
    geom = { geom{:} [3 0.2 0.5] [3 0.2 0.5] [3 1] };

    result = inputgui( geom, uilist, 'pophelp(''pop_biosig'')', ...
                                 'Load data using BIOSIG -- pop_biosig()');
    if length(result) == 0 return; end;
    
    % decode GUI params
    % -----------------
    options = {};
    if ~isempty(result{1}), options = { options{:} 'channels'   eval( [ '[' result{1} ']' ] ) }; end;
    if ~isempty(result{2}), options = { options{:} 'blockrange' eval( [ '[' result{2} ']' ] ) }; end;
    if length(result) > 2
        if ~isempty(result{5}), options = { options{:} 'ref'        eval( [ '[' result{5} ']' ] ) }; end;
        if ~result{3},          options = { options{:} 'rmeventchan' 'off' }; end;
        if  result{4},          options = { options{:} 'blockepoch'  'off'  }; end;
    end;
else
    options = varargin;
end;

% decode imput parameters
% -----------------------
g = finputcheck( options, { 'blockrange'  'integer' [0 Inf]    [];
                            'channels'    'integer' [0 Inf]    [];
                            'ref'         'integer' [0 Inf]    [];
                            'rmeventchan' 'string'  { 'on' 'off' } 'on';
                            'blockepoch'  'string'  { 'on' 'off' } 'off' }, 'pop_biosig');
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
    DAT=sread(dat, newblockrange(2)-newblockrange(1), newblockrange(1));
else 
    DAT=sread(dat, Inf);% this isn't transposed in original!!!!!!!!
    newblockrange    = [];
end
sclose(dat);

if strcmpi(g.blockepoch, 'off')
    dat.NRec = 1;
end;

if ~isempty(newblockrange)
    interval(1) = newblockrange(1) * dat.SampleRate(1) + 1;
    interval(2) = newblockrange(2) * dat.SampleRate(1);
else interval = [];
end
    
EEG = biosig2eeglab(dat, DAT, interval);

if strcmpi(g.rmeventchan, 'on') & strcmpi(dat.TYPE, 'BDF') & isfield(dat, 'BDF')
    if size(EEG.data,1) >= dat.BDF.Status.Channel, 
        disp('Removing event channel...');
        EEG.data(dat.BDF.Status.Channel,:) = []; 
        if ~isempty(EEG.chanlocs)
            EEG.chanlocs(dat.BDF.Status.Channel,:) = [];
        end;
    end;
    EEG.nbchan = size(EEG.data,1);
end;

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
    command = sprintf('EEG = pop_biosig(''%s'');', filename); 
else
    command = sprintf('EEG = pop_biosig(''%s'', %s);', filename, vararg2str(options)); 
end;    
