% eeg_eegrej() - reject porition of continuous data in an EEGLAB 
%                dataset
%
% Usage:
%   >> EEGOUT = eeg_eegrej( EEGIN, regions );
%
% Inputs:
%   INEEG      - input dataset
%   regions    - array of regions to suppress. number x [beg end]  of 
%                regions. 'beg' and 'end' are expressed in term of points
%                in the input dataset. Size of the array is
%                number x 2 of regions.
%
% Outputs:
%   INEEG      - output dataset with updated data, events latencies and 
%                additional boundary events.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 8 August 2002
%
% See also: eeglab(), eegplot(), pop_rejepoch()

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
% Revision 1.22  2004/06/08 20:08:34  arno
% new insertbound
%
% Revision 1.21  2004/06/08 17:09:51  arno
% checking consistency later
%
% Revision 1.20  2004/06/07 18:39:18  arno
% replace indold by indnew
%
% Revision 1.19  2004/06/04 01:30:59  arno
% better msg
%
% Revision 1.18  2004/06/02 17:38:41  arno
% index of index problem
%
% Revision 1.17  2004/06/02 17:26:53  arno
% index to old event
%
% Revision 1.16  2004/06/02 17:21:54  arno
% typo
%
% Revision 1.15  2004/06/02 17:18:00  arno
% insert bound before removing events
%
% Revision 1.14  2004/05/24 20:56:35  arno
% checking event consistency
%
% Revision 1.13  2004/05/14 21:12:51  arno
% same
%
% Revision 1.12  2004/05/14 21:10:10  arno
% update insertbound call
%
% Revision 1.11  2004/04/20 02:08:45  arno
% calling eeg_insertbound
%
% Revision 1.10  2004/03/20 01:48:16  arno
% same
%
% Revision 1.9  2004/03/20 01:46:50  arno
% debug last
%
% Revision 1.8  2004/03/20 01:44:27  arno
% removing bother message of removed events
%
% Revision 1.7  2003/02/27 03:06:47  arno
% debug, handle regions from eegplot
%
% Revision 1.6  2002/11/15 02:14:07  arno
% diabling rejection on
% disk
%
% Revision 1.5  2002/10/22 17:16:53  arno
% debug command line call
%
% Revision 1.4  2002/08/08 14:46:46  arno
% boundary events
%
% Revision 1.3  2002/08/08 01:20:13  arno
% same
%
% Revision 1.2  2002/08/08 01:18:33  arno
% empty regions ?
%
% Revision 1.1  2002/08/08 01:17:31  arno
% Initial revision
%
% Revision 1.12  2002/08/07 22:41:01  arno
% editing text

function [EEG, com] = eeg_eegrej( EEG, regions);

com = '';
if nargin < 2
	help eeg_eegrej;
	return;
end;	
if isempty(regions)
	return;
end;

if isfield(EEG.event, 'latency'), 
	tmpalllatencies = cell2mat( { EEG.event.latency } );
else tmpalllatencies = []; 
end;

% handle regions from eegplot
% ---------------------------
if size(regions,2) > 2, regions = regions(:, 3:4); end;

eeg_options; 
[EEG.data EEG.xmax tmpalllatencies boundevents] = eegrej( EEG.data, ...
												  regions, EEG.xmax-EEG.xmin, tmpalllatencies);
% the string option has been disable since it was causing problems
% ----------------------------------------------------------------
%[EEG.data EEG.xmax tmpalllatencies boundevents] = eegrej( fastif(option_keepdataset, EEG.data, 'EEG.data'), ...
%												  regions, EEG.xmax-EEG.xmin, tmpalllatencies);
EEG.pnts = size(EEG.data,2);
EEG.xmax = EEG.xmax+EEG.xmin;

% change event latencies
% ----------------------
if ~isempty(tmpalllatencies)
	for tmpindex = 1:length(tmpalllatencies)
		EEG.event = setfield(EEG.event, { tmpindex }, 'latency', tmpalllatencies(tmpindex));
	end;

	tmpnanloc = find(isnan(tmpalllatencies));
	EEG.event(tmpnanloc) = [];
    
    if length(boundevents)
        fprintf('eeg_eegrej(): %d boundary (break) events added.\n', length(boundevents));
    end;
    if length(tmpnanloc) > 0
        fprintf('eeg_eegrej(): event latencies recomputed and %d events removed.\n', ...
                length(tmpnanloc));
    end;
end;
EEG.icaact = [];
 
% add boundary events
% -------------------
if ~isempty(boundevents) % boundevent latencies will be recomputed in the function below
    [ EEG.event ] = eeg_insertbound(EEG.event, EEG.pnts, regions);
    EEG = eeg_checkset(EEG, 'eventconsistency');
end;

com = sprintf('%s = eeg_eegrej( %s, %s);', inputname(1), inputname(1), vararg2str({ regions })); 
return;
