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

function res = subsref(obj,s)

    if strcmpi(s(1).type, '.')
        res = builtin('subsref', obj, s);
        return;
    end;
    
    tmpMMO = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });

    subs = s(1).subs;
    finaldim = cellfun('length', subs);
    
    % one dimension input
    % -------------------
    if length(s) > 1 || ~strcmpi(s(1).type, '()')
        error('MMO can only map single array data files');
    end;
    
    % deal with transposed data
    % -------------------------
    if obj.transposed, s = transposeindices(obj, s); end;

    % convert : to real sizes
    % -----------------------
    lastdim = length(subs);
    if isstr(subs{end}) && ndims(obj) > lastdim
        for index = lastdim+1:ndims(obj)
            if index > length(obj.dimensions)
                subs{index} = 1;
            else
                subs{index} = [1:obj.dimensions(index)]; 
            end;
        end;
    end;
    for index = 1:length(subs)
        if isstr(subs{index}) % can only be ":"
            if index > length(obj.dimensions)
                subs{index} = 1;
            else
                subs{index} = [1:obj.dimensions(index)]; 
            end;
        end;
    end;
    finaldim = cellfun(@length, subs);
    finaldim(lastdim) = prod(finaldim(lastdim:end));
    finaldim(lastdim+1:end) = [];

    % non-transposed data
    % -------------------
    res = tmpMMO.data.x(subs{:});
    if length(finaldim) == 1, finaldim(2) = 1; end;
    res = reshape(res, finaldim);
    if obj.transposed
        if finaldim(end) == 1, finaldim(end) = []; end;
        if length(finaldim) <= 2, res = res';
        else
            res = reshape(res, [finaldim(1)*finaldim(2) finaldim(3)])';
            res = reshape(res, [finaldim(3) finaldim(1) finaldim(2)]);
        end;
    end;
    
        
