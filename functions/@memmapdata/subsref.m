% subsref() - index eegdata class
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Nov. 2008

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD
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

function b = subsref(a,s)
    
    if s(1).type == '.'
        b = builtin('subsref', struct(a), s); return;
    end;
        
    subs = s(1).subs;
    finaldim = cellfun('length', subs);
    
    % one dimension input
    % -------------------
    if length(s(1).subs) == 1 
        if isstr(subs{1})
            subs{1} = [1:size(a,1)];
            subs{2} = [1:size(a,2)];
            if ndims(a) == 3, subs{3} = [1:size(a,3)]; end;
            finaldim = prod(size(a));
        end;
        
    % two dimension input
    % -------------------
    elseif length(s(1).subs) == 2 
        if isstr(subs{1}), subs{1} = [1:size(a,1)]; end;
        
        if isstr(subs{2}),
            subs{2} = [1:size(a,2)];
            if ndims(a) == 3, subs{3} = [1:size(a,3)]; end;
        end;
        if length(subs) == 3
             finaldim = [ length(subs{1}) length(subs{2})*length(subs{3}) ];
        else finaldim = [ length(subs{1}) length(subs{2}) ];
        end;
            
    % three dimension input
    % ---------------------
    elseif length(s(1).subs) == 3
        
        if isstr(subs{1}), subs{1} = [1:size(a,1)]; end;
        if isstr(subs{2}), subs{2} = [1:size(a,2)]; end;
        if ndims(a) == 2, 
            subs(3) = []; 
        else
            if isstr(subs{3}), subs{3} = [1:size(a,3)]; end;
        end;
        finaldim = cellfun('length', subs);
     
    end;

    % non-transposed data
    % -------------------
    if ~strcmpi(a.fileformat, 'transposed')
        if length(subs) == 1, b = a.data.data.x(subs{1}); end;
        if length(subs) == 2, b = a.data.data.x(subs{1}, subs{2}); end;
        if length(subs) == 3, b = a.data.data.x(subs{1}, subs{2}, subs{3}); end;
    else
        if ndims(a) == 2
            %if length(s) ==         0, b = transpose(a.data.data.x); return; end;
            if length(s(1).subs) == 1, b = a.data.data.x(s(1).subs{1})'; end;
            if length(s(1).subs) == 2, b = a.data.data.x(s(1).subs{2}, s(1).subs{1})'; end;
            if length(s(1).subs) == 3, b = a.data.data.x(s(1).subs{2}, s(1).subs{1})'; end;
        else
            %if length(s) ==         0, b = permute(a.data.data.x, [3 1 2]); return; end;
            if length(subs) == 1,
                inds1 = mod(subs{1}-1, size(a,1))+1;
                inds2 = mod((subs{1}-inds1)/size(a,1), size(a,2))+1;
                inds3 = ((subs{1}-inds1)/size(a,1)-inds2+1)/size(a,2)+1;
                inds  = (inds1-1)*size(a,2)*size(a,3) + (inds3-1)*size(a,2) + inds2;
                b = a.data.data.x(inds);
            else
                if length(subs) < 2, subs{3} = 1; end;
                
                % repmat if several indices in different dimensions
                % -------------------------------------------------
                len = cellfun('length', subs);
                subs{1} = repmat(reshape(subs{1}, [len(1) 1 1]), [1 len(2) len(3)]);
                subs{2} = repmat(reshape(subs{2}, [1 len(2) 1]), [len(1) 1 len(3)]);
                subs{3} = repmat(reshape(subs{3}, [1 1 len(3)]), [len(1) len(2) 1]);
                
                inds = (subs{1}-1)*a.data.Format{2}(1)*a.data.Format{2}(2) + (subs{3}-1)*a.data.Format{2}(1) + subs{2};
                inds = reshape(inds, [1 prod(size(inds))]);
                b = a.data.data.x(inds);
            end;
        end;
    end;
 
    if length(finaldim) == 1, finaldim(2) = 1; end;
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


