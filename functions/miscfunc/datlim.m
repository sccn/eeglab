% datlim() - return min and max of a matrix
%
% Usage: 
%          >> limits_vector = datlim(data);
%
% Input:
%          data - numeric array
% Outputs:
%          limits_vector = [minval maxval]
%
% Author: Scott Makeig, SCCN/INC/UCSD, May 28, 2005

% Copyright (C) Scott Makeig, SCCN/INC/UCSD, May 28, 2005
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

function [limits_vector] = datlim(data)

if ~isnumeric(data)
   error('data must be a numeric array')
   return
end

limits_vector = [ min(data(:)) max(data(:)) ]; % thanks to Arno Delorme

% minval = squeeze(min(data)); maxval = squeeze(max(data));
% while numel(minval) > 1
%    minval = squeeze(min(minval)); maxval = squeeze(max(maxval));
% end
% limits_vector = [minval maxval];
   

