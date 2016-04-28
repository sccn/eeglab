% eegthresh() -  reject trials with out-of-bounds channel values within a
%                specified epoch time range.
%
% Usage:
%   >> [Iin, Iout, newsignal, elec] = eegthresh( signal, frames, ...
%                elecs, negthresh, posthresh, timerange, starttime,endtime);
%
% Required inputs:
%   signal     - 2-D data matrix [channels, frames*sweeps],  
%                or 3-D data matrix [channels, frames, sweeps]
%   frames     - number of points per epoch
%   elecs      - [int vector] electrode indices to reject on
%   negthresh  - minimum rejection threshold(s) in uV. This can be an array 
%                of values for individual electrodes. If fewer values than  
%                electrodes, the last value is used for the remaining 
%                electrodes.
%   posthresh  - maximum rejection threshold(s) in uV (same syntax as for
%                negthresh)
%   timerange  - [mintime maxtime] time range limits of the signal
%   starttime  - Starting limit (in seconds or Hz) of range to perform 
%                rejection (same syntax  as negthresh)
%   endtime    - Ending limit (in seconds or Hz) of range to perform 
%                rejection (same syntax  as negthresh)
%
% Outputs:
%   Iin        - Indexes of epochs accepted
%   Iout       - Indexes of epochs rejected
%   newsignal  - input data after epoch rejection
%   elec       - electrode that triggered the rejections (array of 0s 
%                and 1s with the same number of columns as Iout
%                and number of rows = number of electrodes).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_eegthresh(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [Iin, Iout, newsignal, elec] = eegthresh( signal, pnts, electrodes, negthresh, posthresh, timerange, starttime, endtime);

if nargin < 7
	help eegthresh;
	return;
end;

if starttime < timerange(1)
	disp('eegthresh: starting point out of range, adjusted');
    starttime = timerange(1);
end;	
if endtime > timerange(2)
	disp('eegthresh: ending point out of range, adjusted');
    endtime = timerange(2);
end;	

% complete thresholds values if necessary
%----------------------------------------
if size(posthresh,2) < size(electrodes,2)
	posthresh = [ posthresh posthresh(end)*ones(1,size(electrodes,2)-size(posthresh,2))];
end;	
if size(negthresh,2) < size(electrodes,2)
	negthresh = [ negthresh negthresh(end)*ones(1,size(electrodes,2)-size(negthresh,2))];
end;
	
% complete timeselect values if necessary
%----------------------------------------
if size(starttime,2) < size(electrodes,2)
	starttime = [ starttime starttime(end)*ones(1,size(electrodes,2)-size(starttime,2))];
end;	
if size(endtime,2) < size(electrodes,2)
	endtime = [ endtime endtime(end)*ones(1,size(electrodes,2)-size(endtime,2))];
end;	

% find the maximum for each trial
%--------------------------------
sweeps = size(signal(1,:),2)/pnts;
signal = reshape(signal(electrodes,:), size(electrodes(:),1), pnts, sweeps);

% reject the selected trials
%---------------------------
elec = zeros(size(electrodes(:),1), sweeps);
allelec = zeros(1, sweeps);
for indexe = 1:size(electrodes(:),1)
	% transform the time range
	% ------------------------
	framelowlimit  = max(1,floor((starttime(indexe)-timerange(1))/(timerange(2)-timerange(1))*(pnts-1))+1);	
	framehighlimit = floor((endtime(indexe)  -timerange(1))/(timerange(2)-timerange(1))*(pnts-1))+1;

	% remove outliers
	% ---------------
	sigtmp = squeeze(signal(indexe,framelowlimit:framehighlimit,:));
    if size(signal,3) == 1, sigtmp = sigtmp'; end;
	sigmax = max(sigtmp, [], 1);
	sigmin = min(sigtmp, [], 1);
	elec(indexe,:) = ( sigmin < negthresh(indexe) ) | ( sigmax > posthresh(indexe) );
	allelec = allelec | elec(indexe,:);
end;
Iout = find( allelec == 1 );
Iin  = find( allelec == 0 );
elec = elec(:,Iout);

% reject the selected trials
%---------------------------
newsignal = reshape(signal, size(signal,1), pnts, sweeps);
if ~isempty(Iin)
	newsignal = newsignal(:,:,Iin);
	if ndims(signal) == 2
		newsignal = newsignal(:,:);
	end;
else
	newsignal = [];
end;
return;
