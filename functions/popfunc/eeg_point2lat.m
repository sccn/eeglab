% eeg_point2lat() - convert latency in data points to latency in ms relative
%                   to the time locking. Used in eeglab().
% Usage:
%       >> [newlat ] = eeg_point2lat( lat_array, epoch_array,...
%                                 srate, timelimits, timeunit);
% Inputs:
%   lat_array   - latency array in data points assuming concatenated
%                 data epochs (see eeglab() event structure)
%   epoch_array - epoch number corresponding to each latency value
%   srate       - data sampling rate in Hz
%   timelimits  - [min max] timelimits in 'timeunit' units (see below)
%   timeunit    - time unit in second. Default is 1 = seconds.
%
% Outputs:
%   newlat      - converted latency values (in 'timeunit' units) for each epoch
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2 Mai 2002
%
% See also: eeg_point2lat(), eeglab(), pop_editieventvals(), pop_loaddat()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2 Mai 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.4  2003/01/01 18:23:57  scott
% header edit & unused output noted
%
% Revision 1.3  2002/08/08 22:03:33  arno
% update
%
% Revision 1.2  2002/08/08 14:55:55  arno
% rounding very low latencies
%
% Revision 1.1  2002/05/03 01:04:42  arno
% Initial revision
%

function newlat = eeg_point2lat( lat_array, epoch_array, srate, timewin, timeunit);

if nargin <4
    help eeg_point2lat;
    return;
end;    
if nargin <5
	timeunit = 1;
end;

if length(lat_array) ~= length(epoch_array)
	if length(epoch_array)~= 1
		disp('eeg_point2lat: latency and epoch arrays must have the same length'); return;
	else
		epoch_array = ones(1,length(lat_array))*epoch_array;
	end;
end;
if length(timewin) ~= 2
    disp('eeg_point2lat: timelimits array must have length 2'); return;
end;
if iscell(epoch_array)
	epoch_array = cell2mat(epoch_array);
end;
if iscell(lat_array)
	lat_array = cell2mat(lat_array);
end

timewin = timewin*timeunit;

pnts = (timewin(2)-timewin(1))*srate+1;
newlat  = ((lat_array - (epoch_array-1)*pnts-1)/srate+timewin(1))/timeunit;
newlat = round(newlat*1E9)*1E-9;
