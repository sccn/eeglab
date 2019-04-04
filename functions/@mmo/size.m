% size() - size of memory mapped underlying array
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Nov. 2008

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD
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

function [s varargout] = size(obj,dim)
    
    if obj.transposed
        if length(obj.dimensions) ~= 2 && length(obj.dimensions) ~= 3
            error('Transposed array must be 2 or 3 dims');
        end
        if length(obj.dimensions) == 2 tmpdimensions = [obj.dimensions(2) obj.dimensions(1)];
        else                           tmpdimensions = [obj.dimensions(3) obj.dimensions(1:end-1)];
        end
    else
        tmpdimensions = obj.dimensions;
    end
    
    s = tmpdimensions;
    
    if nargin > 1
        if dim >length(s)
            s = 1;
        else
            s = s(dim);
        end
    else
        if nargout > 1
            s = [s 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
            alls = s;
            s = s(1);
            for index = 1:max(nargout,1)-1
                varargout{index} = alls(index+1);
            end
        end
    end
    
    
