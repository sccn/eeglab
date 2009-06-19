% pop_readbdf() - Read Biosemi 24-bit BDF file 
%
% Usage:
%   >> EEG = pop_readbdf;             % an interactive window pops up
%   >> EEG = pop_readbdf( filename ); % no pop-up window 
%   >> EEG = pop_readbdf( filename, range, eventchans, ref );
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
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 13 March 2002
%
% See also: openbdf(), readbdf()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 13 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% programmed from pop_readedf() version 1.15

% $Log: not supported by cvs2svn $
% Revision 1.24  2004/05/24 21:44:32  arno
% same
%
% Revision 1.23  2004/05/24 21:43:47  arno
% show number of channels
%
% Revision 1.22  2004/05/15 22:56:50  scott
% edit text
%
% Revision 1.21  2004/05/08 16:06:01  scott
% edited text box, help message
%
% Revision 1.20  2004/03/20 00:51:25  arno
% message
%
% Revision 1.19  2004/03/19 00:50:25  arno
% header msg
%
% Revision 1.18  2004/03/18 02:09:13  arno
% change channel event default
%
% Revision 1.17  2004/01/29 02:05:03  arno
% gui edit
%
% Revision 1.16  2004/01/29 02:02:53  arno
% putting block number in GUI
%
% Revision 1.15  2003/12/18 00:28:46  arno
% reading channel label
%
% Revision 1.14  2003/12/18 00:10:10  arno
% same
%
% Revision 1.13  2003/12/18 00:08:27  arno
% data bug
%
% Revision 1.12  2003/12/18 00:03:37  arno
% debug
% last
%
% Revision 1.11  2003/12/18 00:02:00  arno
% remove event channels
%
% Revision 1.10  2003/12/17 19:01:16  arno
% do not remove reference channels
%
% Revision 1.9  2003/12/17 17:48:40  arno
% adding entry for reference
%
% Revision 1.8  2003/12/11 20:07:59  arno
% eventindices -> eventchans
% .,
%
% Revision 1.7  2003/12/11 20:01:48  arno
% debug history
%
% Revision 1.6  2003/12/11 19:23:50  arno
% header and now import events automatically
%
% Revision 1.5  2003/11/04 01:22:40  arno
% removing warnings
%
% Revision 1.4  2003/10/13 17:42:44  arno
% fix blockrage
%
% Revision 1.3  2003/10/10 17:11:12  arno
% fixing default arg
%
% Revision 1.2  2003/08/29 18:44:16  arno
% adding gui for block range
%
% Revision 1.1  2003/06/05 15:38:44  arno
% Initial revision
%

function [EEG, command] = pop_readbdf(filename, blockrange, eventchans, ref); 
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
                     'If Biosemi data, reference chan(s) number(s)' };
    inistr       = { '' int2str(dat.Head.NS) '' };
    result       = inputdlg2( promptstr, 'Import BDF file -- pop_readbdf()', 1,  inistr, 'pop_readbdf');
    if length(result) == 0 return; end;
    blockrange   = eval( [ '[' result{1} ']' ] );
    eventchans   = eval( [ '[' result{2} ']' ] );
    ref          = eval( [ '[' result{3} ']' ] );
end;

% load datas
% ----------
EEG = eeg_emptyset;
fprintf('Reading BDF data in 24-bit format...\n');
dat = openbdf(filename);
%dat.Head = sopen(filename);
dat.Head.NRec = 100;
if isempty(dat.Head.NRec)
    dsafsd
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
if ~isempty(eventchans)
    if any(eventchans > EEG.nbchan) | any(eventchans < 1)
        error('Event channel index or indices out of range');
    end;
    disp('Applying 8-bit masks to event channels (only for BIOSEMI files).');
    EEG.data(eventchans,:) = bitand(uint16(EEG.data(eventchans,:)), 255);
    EEG = pop_chanevent(EEG, eventchans, 'edge', 'leading', 'delchan', 'on');
end;

% rerefencing
% -----------
disp('Re-referencing...');
if ~isempty(ref)
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
