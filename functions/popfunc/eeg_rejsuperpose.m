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

if typerej == 0 | Rothertype
	if Rmanual
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejmanual); % see bottom for the
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejmanualE); % function rejarray
	end;
	if Rthres
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejthresh);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejthreshE);
	end;
	if Rfreq
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejfreq);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejfreqE);
	end;
	if Rconst
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejconst);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejconstE);
	end;
	if Rent
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejjp);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejjpE);
	end;
	if Rkurt
		rejglobal  = rejarray( rejglobal,  EEG.reject.rejkurt);
		rejglobalE = rejarray( rejglobalE, EEG.reject.rejkurtE);
	end;
end;

% ---------------

if typerej == 1 | Rothertype
	if Rmanual
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejmanual);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejmanualE);
	end;
	if Rthres
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejthresh);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejthreshE);
	end;
	if Rfreq
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejfreq);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejfreqE);
	end;
	if Rconst
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejconst);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejconstE);
	end;
	if Rent
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejjp);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejjpE);
	end;
	if Rkurt
		rejglobal  = rejarray( rejglobal,  EEG.reject.icarejkurt);
		rejglobalE = rejarray( rejglobalE, EEG.reject.icarejkurtE);
	end;
end;

EEG.reject.rejglobal = rejglobal;
EEG.reject.rejglobalE = rejglobalE;

com =sprintf('%s = eeg_rejsuperpose( %s, %d, %d, %d, %d, %d, %d, %d, %d);', ...
	  inputname(1), inputname(1), ~typerej, Rmanual, Rthres, Rconst, Rent, Rkurt, Rfreq, Rothertype);

return;

% subfunction rejecting an array ------ 
function dest = rejarray( dest, ori)
    if isempty(dest)
        dest = ori;
    elseif ~isempty(ori)
		dest = dest | ori;
	end;
return;
