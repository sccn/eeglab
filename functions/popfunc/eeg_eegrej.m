% eeg_eegrej() - reject porition of continuous data in an EEGLAB 
%                dataset
%
% Usage:
%   >> EEGOUT = eeg_eegrej( EEGIN, regions );
%
% Inputs:
%   INEEG      - input dataset
%   regions    - array of regions to suppress. [beg end] x number of 
%                regions. 'beg' and 'end' are expressed in term of points
%                in the input dataset. Size of the array is
%                2xnumber of regions.
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

eeg_options; 
[EEG.data EEG.xmax tmpalllatencies boundevents] = eegrej( fastif(option_keepdataset, EEG.data, 'EEG.data'), ...
												  regions(:,3:4), EEG.xmax-EEG.xmin, tmpalllatencies);
EEG.pnts = size(EEG.data,2);
EEG.xmax = EEG.xmax+EEG.xmin;

% change event latencies
% ----------------------
if ~isempty(tmpalllatencies)
	tmpnanloc = find(~isnan(tmpalllatencies));
	EEG.event = EEG.event(tmpnanloc);
	fprintf('eeg_eegrej(): event latencies recomputed and %d (of %d) events removed.\n', ...
			length(tmpalllatencies)-length(EEG.event), length(tmpalllatencies));
	tmpalllatencies = tmpalllatencies(tmpnanloc);
	for tmpindex = 1:length(EEG.event)
		EEG.event = setfield(EEG.event, { tmpindex }, 'latency', tmpalllatencies(tmpindex));
	end;
end;

% add boundary events
% -------------------
if ~isempty(boundevents)
	fprintf('eeg_eegrej(): boundary events added.\n');
	for tmpindex = 1:length(EEG.event)
		EEG.event(end+1).type  = 'boundary';
		EEG.event(end).latency = boundevents(tmpindex);
	end;
	EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
	EEG = eeg_checkset(EEG, 'eventconsistency');
end;
EEG.icaact = [];

com = sprintf('%s = eeg_eegrej( %s, %s);', inputname(1), inputname(1), vararg2str({ regions(:,1:4) })); 
return;
