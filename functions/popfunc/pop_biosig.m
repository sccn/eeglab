% pop_biosig() - import data files into EEGLAB using BIOSIG toolbox
%
% Usage:
%   >> OUTEEG = pop_biosig; % pop up window
%   >> OUTEEG = pop_biosig( filename, channels, type);
%
% Inputs:
%   filename - [string] file name
%   channels - [integer array] list of channel indices
%   type     - [string] file type. See sload() function. Not functional
%              because sload() does not support it anymore.
%
% Outputs:
%   OUTEEG   - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Oct. 29, 2003
%
% Note: BIOSIG toolbox must be installed. Download BIOSIG at 
%       http://sourceforge.net/project/showfiles.php?group_id=7072
%       (please let us know at eeglab@sccn.ucsd.edu if this link is outdated).
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
% Revision 1.3  2005/10/27 05:24:55  arno
% filename
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

function [EEG, command] = pop_biosig(filename, channels, type); 
EEG = [];
command = '';

if nargin < 1
    if ~exist('inputgui')
        error('requires the EEGLAB toolbox for pop up window');
    end;
    
	% ask user
	[filename, filepath] = uigetfile('*.*', 'Choose a DAQ file -- pop_biosig'); 
    drawnow;
	if filename == 0 return; end;
	filename = fullfile(filepath,filename);
    %alltypes = { 'autodetect' 'DAQ' 'MAT' 'DAT' 'EBS' 'RAW' 'RDT' 'XLS' 'DA_' 'RG64' 'SIG' };

    uilist   = { { 'style' 'text' 'String' 'Channel list (defaut all):' } ...
                 { 'style' 'edit' 'string' '' } };
    %             { 'style' 'text' 'string' 'Force data type' } ...
    %             { 'style' 'list' 'string'  strvcat(alltypes{:}) } };
    
    results = inputgui( { [1 1] }, uilist, 'pophelp(''pop_biosig'')', ...
                                 'Load data using BIOSIG -- pop_biosig()');
    if length(results) == 0 return; end;
    
    if ~isempty(results{1})
        channels = eval( [ '[ ' results{1} ' ]' ]);
    end;
    %if results{2} ~= 1
    %    type = alltypes{results{2}};
    %end;
end;

% loading data
% ------------
if exist('channels') ~= 1
    channels = 0;
end;
disp('Importing data...');
[signal,H] = sload(filename,channels);
        
% decoding data
% -------------
try, EEG = eeg_emptyset;
catch, end;
[signal, H]  = sload(filename);
EEG.comments = [ 'Original file: ' filename ];
EEG.srate    = H.SampleRate(1);
EEG.data     = signal';
EEG.nbchan   = size(EEG.data,1);
EEG.trials   = H.NRec;
EEG.pnts     = size(EEG.data,2)/H.NRec;
if isfield(H, 'Label') & ~isempty(H.Label)
    EEG.chanlocs        = struct('labels', cellstr(char(H.Label)));
end;
if isfield(H, 'EVENT')    
    EEG.event = biosig2eeglabevent(H.EVENT);
else 
    disp('Warning: no event found. Events might be embeded in a data channel.');
    disp('         To extract events, use menu File > Import Event Info > From data channel');
end;

if exist('type') ~= 1
    command = sprintf('EEG = pop_biosig(''%s'', %s);', filename, vararg2str({ channels }));
else
    command = sprintf('EEG = pop_biosig(''%s'', %s);',filename, vararg2str({ channels type }));
end;
return;
