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
   

