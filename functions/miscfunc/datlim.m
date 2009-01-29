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
   

