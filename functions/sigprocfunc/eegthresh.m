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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [Iin, Iout, newsignal, elec] = eegthresh( signal, pnts, electrodes, negthresh, posthresh, timerange, starttime, endtime);

if nargin < 7
	help eegthresh;
	return;
end

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
end
	
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
    if size(signal,3) == 1, sigtmp = sigtmp'; end
	sigmax = max(sigtmp, [], 1);
	sigmin = min(sigtmp, [], 1);
	elec(indexe,:) = ( sigmin < negthresh(indexe) ) | ( sigmax > posthresh(indexe) );
	allelec = allelec | elec(indexe,:);
end
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
	end
else
	newsignal = [];
end
return;
