% eeg_rejsuperpose() - superpose rejections of a EEG dataset.
%
% Usage:
%   >> EEGOUT = eeg_rejsuperpose( EEGIN, typerej, Rmanual, Rthres, ...
%                    Rconst, Rent, Rkurt, Rfreq, Rothertype);
%
% Inputs:
%   EEGIN      - input dataset
%   typerej    - type of rejection (1=raw data; 0=ica).
%   Rmanual    - include manual rejection (0|1). 
%   Rthres     - include threshold rejection (0|1). 
%   Rconst     - include rejection of constant activity (0|1). 
%   Rent       - include entropy rejection (0|1). 
%   Rkurt      - include kurtosis rejection (0|1). 
%   Rfreq      - include frequcy based rejection (0|1). 
%   Rothertype - include manual rejection (0|1). 
%
% Outputs:
%   EEGOUT     - with rejglobal and rejglobalE fields updated
%
% See also: eeglab()

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

function [EEG, com] = eeg_rejsuperpose( EEG, typerej, Rmanual, Rthres, Rconst, ...
                              Rent, Rkurt, Rfreq, Rothertype);

if nargin < 9
	help eeg_rejsuperpose;
	return;
end;	

typerej = ~typerej;
rejglobal  = zeros( 1, EEG.trials);
if typerej == 0
    rejglobalE = zeros( EEG.nbchan, EEG.trials);
else
    rejglobalE = zeros( size(EEG.icaweights,1), EEG.trials);
end

if typerej == 0 || Rothertype
	if Rmanual
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejmanual); % see bottom for the
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejmanualE); % function rejarray
	end
	if Rthres
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejthresh);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejthreshE);
	end
	if Rfreq
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejfreq);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejfreqE);
	end
	if Rconst
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejconst);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejconstE);
	end
	if Rent
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejjp);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejjpE);
	end
	if Rkurt
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejkurt);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejkurtE);
	end
end

% ---------------

if typerej == 1 || Rothertype
	if Rmanual
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejmanual);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejmanualE);
	end
	if Rthres
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejthresh);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejthreshE);
	end
	if Rfreq
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejfreq);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejfreqE);
	end
	if Rconst
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejconst);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejconstE);
	end
	if Rent
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejjp);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejjpE);
	end
	if Rkurt
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejkurt);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejkurtE);
	end
end

EEG.reject.rejglobal = rejglobal;
EEG.reject.rejglobalE = rejglobalE;

com =sprintf('EEG = eeg_rejsuperpose( EEG, %d, %d, %d, %d, %d, %d, %d, %d);', ...
	   ~typerej, Rmanual, Rthres, Rconst, Rent, Rkurt, Rfreq, Rothertype);

return;

% subfunction rejecting an array ------ 
function dest = rejarray( dest, ori)
    if isempty(dest)
        dest = ori;
    elseif ~isempty(ori)
		dest = dest | ori;
	end
return;
