% pop_read_erpss() - interactively import an uncompressed ERPSS-format 
%                    data file (.RAW or .RDF) using read_erpss()
% Usage:
%   >> OUTEEG = pop_read_erpss( filename, srate );
%
% Inputs:
%   filename       - file name (with extension)
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 23 January 2003
%
% See also: eeglab(), read_erpss()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.7  2003/04/10 18:00:15  arno
% file file
% file filter
%
% Revision 1.6  2003/01/28 18:15:53  arno
% adding filename and filepath
%
% Revision 1.5  2003/01/28 16:49:09  arno
% debugging read error
%
% Revision 1.4  2003/01/24 17:49:39  arno
% now reading sampling rate
%
% Revision 1.3  2003/01/24 04:13:00  scott
% header edit -sm
%
% Revision 1.2  2003/01/24 01:34:44  arno
% adding setname
%
% Revision 1.1  2003/01/24 01:29:57  arno
% Initial revision
%

function [EEG, command] = pop_read_erpss(filename); 
command = '';

if nargin < 1
	% ask user
	[filenametmp, filepath] = uigetfile('*.rdf;*.RDF;*.raw;*.RAW', 'Choose a ERPSS file -- pop_read_erpss'); 
    drawnow;
	if filenametmp == 0 return; end;
	filename = [filepath '/' filenametmp];
end;

% read ERPSS format
EEG = eeg_emptyset;
fprintf('pop_read_erpss: importing ERPSS file...\n');
[EEG.data,events,header] = read_erpss(filename);
EEG.nbchan = size(EEG.data,1);
EEG.srate = header.srate;
EEG.setname = 'ERPSS data';
if exist('filepath') == 1
    EEG.filepath = filepath;
    EEG.filename = filenametmp;
else
    EEG.filename = filename;
end;
EEG.event = struct( 'type', { events.event_code }, 'latency', {events.sample_offset});
EEG = eeg_checkset(EEG);

command = sprintf('EEG = pop_read_erpss(''%s'');',filename); 
return;
