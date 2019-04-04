% subsref() - index eegdata class
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

function res = subsref(obj,s)

    if strcmpi(s(1).type, '.')
        res = builtin('subsref', obj, s);
        return;
    end
    
    tmpMMO = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });

    subs = s(1).subs;
    finaldim = cellfun('length', subs);
    
    % one dimension input
    % -------------------
    if length(s) > 1 || ~strcmpi(s(1).type, '()')
        error('MMO can only map single array data files');
    end
    
    % deal with transposed data
    % -------------------------
    if obj.transposed, s = transposeindices(obj, s); end

    % convert : to real sizes
    % -----------------------
    lastdim = length(subs);
    if ischar(subs{end}) && ndims(obj) > lastdim
        for index = lastdim+1:ndims(obj)
            if index > length(obj.dimensions)
                subs{index} = 1;
            else
                subs{index} = [1:obj.dimensions(index)]; 
            end
        end
    end
    for index = 1:length(subs)
        if ischar(subs{index}) % can only be ":"
            if index > length(obj.dimensions)
                subs{index} = 1;
            else
                subs{index} = [1:obj.dimensions(index)]; 
            end
        end
    end
    finaldim = cellfun(@length, subs);
    finaldim(lastdim) = prod(finaldim(lastdim:end));
    finaldim(lastdim+1:end) = [];

    % non-transposed data
    % -------------------
    res = tmpMMO.data.x(subs{:});
    if length(finaldim) == 1, finaldim(2) = 1; end
    res = reshape(res, finaldim);
    if obj.transposed
        if finaldim(end) == 1, finaldim(end) = []; end
        if length(finaldim) <= 2, res = res';
        else
            res = reshape(res, [finaldim(1)*finaldim(2) finaldim(3)])';
            res = reshape(res, [finaldim(3) finaldim(1) finaldim(2)]);
        end
    end
    
        
