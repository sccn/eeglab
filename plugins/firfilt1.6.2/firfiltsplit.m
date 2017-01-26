% firfiltsplit() - Split data at discontinuities and forward to dc padded
%                  filter function
%
% Usage:
%   >> EEG = firfiltsplit(EEG, b);
%
% Inputs:
%   EEG           - EEGLAB EEG structure
%   b             - vector of filter coefficients
%   causal        - scalar boolean perform causal filtering {default 0}
%
% Outputs:
%   EEG           - EEGLAB EEG structure
%
% Note:
%   This function is (in combination with firfiltdcpadded) just a
%   non-memory optimized version of the firfilt function allowing causal
%   filtering. Will possibly replace firfilt in the future.
%
% Author: Andreas Widmann, University of Leipzig, 2013
%
% See also:
%   firfiltdcpadded, findboundaries

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2013 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function EEG = firfiltsplit(EEG, b, causal)

if nargin < 3 || isempty(causal)
    causal = 0;
end
if nargin < 2
    error('Not enough input arguments.');
end

% Find data discontinuities and reshape epoched data
if EEG.trials > 1 % Epoched data
    EEG.data = reshape(EEG.data, [EEG.nbchan EEG.pnts * EEG.trials]);
    dcArray = 1 : EEG.pnts : EEG.pnts * (EEG.trials + 1);
else % Continuous data
    dcArray = [findboundaries(EEG.event) EEG.pnts + 1];
end

% Loop over continuous segments
for iDc = 1:(length(dcArray) - 1)

    % Filter segment
    EEG.data(:, dcArray(iDc):dcArray(iDc + 1) - 1) = firfiltdcpadded(b, EEG.data(:, dcArray(iDc):dcArray(iDc + 1) - 1)', causal)';

end

% Reshape epoched data
if EEG.trials > 1
    EEG.data = reshape(EEG.data, [EEG.nbchan EEG.pnts EEG.trials]);
end

end
