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

function b = subsref(a,s)
    
    if s(1).type == '.'
        b = builtin('subsref', struct(a), s); return;
    end
        
    subs = s(1).subs;
    finaldim = cellfun('length', subs);
    
    % one dimension input
    % -------------------
    if length(s(1).subs) == 1 
        if ischar(subs{1})
            subs{1} = [1:size(a,1)];
            subs{2} = [1:size(a,2)];
            if ndims(a) == 3, subs{3} = [1:size(a,3)]; end
            finaldim = prod(size(a));
        end
        
    % two dimension input
    % -------------------
    elseif length(s(1).subs) == 2 
        if ischar(subs{1}), subs{1} = [1:size(a,1)]; end
        
        if ischar(subs{2}),
            subs{2} = [1:size(a,2)];
            if ndims(a) == 3, subs{3} = [1:size(a,3)]; end
        end
        if length(subs) == 3
             finaldim = [ length(subs{1}) length(subs{2})*length(subs{3}) ];
        else finaldim = [ length(subs{1}) length(subs{2}) ];
        end
            
    % three dimension input
    % ---------------------
    elseif length(s(1).subs) == 3
        
        if ischar(subs{1}), subs{1} = [1:size(a,1)]; end
        if ischar(subs{2}), subs{2} = [1:size(a,2)]; end
        if ndims(a) == 2, 
            subs(3) = []; 
        else
            if ischar(subs{3}), subs{3} = [1:size(a,3)]; end
        end
        finaldim = cellfun('length', subs);
     
    end

    % non-transposed data
    % -------------------
    if ~strcmpi(a.fileformat, 'transposed')
        if length(subs) == 1, b = a.data.data.x(subs{1}); end
        if length(subs) == 2, b = a.data.data.x(subs{1}, subs{2}); end
        if length(subs) == 3, b = a.data.data.x(subs{1}, subs{2}, subs{3}); end
    else
        if ndims(a) == 2
            %if length(s) ==         0, b = transpose(a.data.data.x); return; end
            if length(s(1).subs) == 1, b = a.data.data.x(s(1).subs{1})'; end
            if length(s(1).subs) == 2, b = a.data.data.x(s(1).subs{2}, s(1).subs{1})'; end
            if length(s(1).subs) == 3, b = a.data.data.x(s(1).subs{2}, s(1).subs{1})'; end
        else
            %if length(s) ==         0, b = permute(a.data.data.x, [3 1 2]); return; end
            if length(subs) == 1,
                inds1 = mod(subs{1}-1, size(a,1))+1;
                inds2 = mod((subs{1}-inds1)/size(a,1), size(a,2))+1;
                inds3 = ((subs{1}-inds1)/size(a,1)-inds2+1)/size(a,2)+1;
                inds  = (inds1-1)*size(a,2)*size(a,3) + (inds3-1)*size(a,2) + inds2;
                b = a.data.data.x(inds);
            else
                if length(subs) < 2, subs{3} = 1; end
                
                % repmat if several indices in different dimensions
                % -------------------------------------------------
                len = cellfun('length', subs);
                subs{1} = repmat(reshape(subs{1}, [len(1) 1 1]), [1 len(2) len(3)]);
                subs{2} = repmat(reshape(subs{2}, [1 len(2) 1]), [len(1) 1 len(3)]);
                subs{3} = repmat(reshape(subs{3}, [1 1 len(3)]), [len(1) len(2) 1]);
                
                inds = (subs{1}-1)*a.data.Format{2}(1)*a.data.Format{2}(2) + (subs{3}-1)*a.data.Format{2}(1) + subs{2};
                inds = reshape(inds, [1 prod(size(inds))]);
                b = a.data.data.x(inds);
            end
        end
    end
 
    if length(finaldim) == 1, finaldim(2) = 1; end
    b = reshape(b, finaldim);

% 2 dims
%inds1 = mod(myinds-1, size(a,1))+1;
%inds2 = (myinds-inds1)/size(a,1)+1;
%inds  = (inds2-1)*size(a,1) + inds1;

% 3 dims
%inds1 = mod(myinds-1, size(a,1))+1;
%inds2 = mod((myinds-inds1)/size(a,1), size(a,2))+1;
%inds3 = ((myinds-inds1)/size(a,1)-inds2)/size(a,2)+1;
%inds  = (inds3-1)*size(a,1)*size(a,2) + inds2*size(a,1) + inds1;


