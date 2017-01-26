% firfiltdcpadded() - Pad data with DC constant and filter
%
% Usage:
%   >> data = firfiltdcpadded(data, b, causal);
%
% Inputs:
%   data      - raw data
%   b         - vector of filter coefficients
%   causal    - boolean perform causal filtering {default 0}
%
% Outputs:
%   data      - smoothed data
%
% Note:
%   firfiltdcpadded always operates (pads, filters) along first dimension.
%   Not memory optimized.
%
% Author: Andreas Widmann, University of Leipzig, 2013
%
% See also:
%   firfiltsplit

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

function [ data ] = firfiltdcpadded(b, data, causal)

% Defaults
if nargin < 3 || isempty(causal)
    causal = 0;
end

% Check arguments
if nargin < 2
    error('Not enough input arguments.');
end

% Filter's group delay
if mod(length(b), 2) ~= 1
    error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;
b = double(b); % Filter with double precision

% Pad data with DC constant
if causal
    startPad = repmat(data(1, :), [2 * groupDelay 1]);
    endPad = [];
else
    startPad = repmat(data(1, :), [groupDelay 1]);
    endPad = repmat(data(end, :), [groupDelay 1]);
end

% Filter data
data = filter(b, 1, double([startPad; data; endPad])); % Pad and filter with double precision

% Remove padded data
data = data(2 * groupDelay + 1:end, :);
 
end
