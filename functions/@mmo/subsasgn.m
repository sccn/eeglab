% subsasgn() - define index assignment for eegdata objects
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

function obj = subsasgn(obj,ss,val)

% check empty assignment
% -----------------------
for index = 1:length(ss(1).subs)
    if isempty(ss(1).subs{index}), return; end
end

% remove useless ":"
% ------------------
while length(obj.dimensions) < length(ss(1).subs)
    if isequal(ss(1).subs{end}, ':')
        ss(1).subs(end) = [];
    else
        break;
    end
end

% deal with transposed data
% -------------------------
if obj.transposed, ss = transposeindices(obj, ss); end

% check stack and local variables
% ---------------------
ncopies = checkworkspace(obj);
if ncopies < 2
    s = evalin('caller', 'whos');
    for index = 1:length(s)
        if strcmpi(s(index).class, 'struct') || strcmpi(s(index).class, 'cell')
            tmpVar = evalin('caller', s(index).name);
            ncopies = ncopies + checkcopies_local(obj, tmpVar);
        elseif strcmpi(s(index).class, 'mmo')
            if s(index).persistent || s(index).global
                disp('Warning: mmo objects should not be made persistent or global. Behavior is unpredictable.');
            end
            tmpVar = evalin('caller', s(index).name);
            if isequal(tmpVar, obj), ncopies = ncopies + 1; end
            if ncopies > 1, break; end
        end
    end
end

% removing some entries
% ---------------------
if isempty(val)
    newdim1 = obj.dimensions;
    newdim2 = obj.dimensions;
    
    % create new array of new size
    % ----------------------------
    tmpMMO = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });
    newFileName = mmo.getnewfilename;
    
    % find non singleton dimension
    % ----------------------------
    nonSingleton = [];
    ss2 = ss;
    for index = 1:length(ss(1).subs)
        if ~ischar(ss(1).subs{index}) % can only be ":"
            nonSingleton(end+1) = index;
            if islogical( ss(1).subs{index}),  ss(1).subs{index} = find( ss(1).subs{index}); end
            ss2(1).subs{index} = setdiff_bc([1:newdim1(index)], ss(1).subs{index}); % invert selection
        end
    end
    if length(nonSingleton) > 1, error('A null assignment can have only one non-colon index'); end
    if isempty(nonSingleton), obj = []; return; end
    
    % compute new final dimension and copy data
    % -----------------------------------------
    if length(ss(1).subs) == 1
        fid = fopen(newFileName, 'w');
        newdim2 = [ prod(newdim2)-length(ss(1).subs{1}) ];
        if ~(newdim1(1) > 1 && all(newdim1(2:end) == 1)), newdim2 = [1 newdim2]; 
        else                                              newdim2 = [newdim2 1]; end
        newindices = setdiff_bc([1:prod(newdim1)], ss(1).subs{1});
        for index = newindices
            fwrite(fid, tmpMMO.Data.x(index), 'float');
        end
        fclose(fid);
        fprintf('Warning: memory mapped object writing might not be up to date in cache on network drive');
    else
        if length(ss(1).subs) < length(newdim2)
            newdim2(length(ss(1).subs)) = prod(newdim2(length(ss(1).subs):end));
            newdim2(length(ss(1).subs)+1:end) = [];
            if nonSingleton == length(ss(1).subs)
                ss2(1).subs{end} = setdiff_bc([1:newdim2(end)], ss(1).subs{end});
            end
        end
        newdim2(nonSingleton) = newdim2(nonSingleton)-length(ss(1).subs{nonSingleton});
        if ischar(ss2.subs{end})
            ss2.subs{end} = [1:prod(newdim1(length(ss2.subs):end))];
        end
        ss3 = ss2;
        fid = fopen(newFileName, 'w');
        
        % copy large blocks
        alllen = cellfun(@length, ss2.subs);
        inc = ceil(250000/prod(alllen(1:end-1))); % 1Mb block
        for index = 1:inc:length(ss2.subs{end})
            ss3.subs{end} = ss2.subs{end}(index:min(index+inc, length(ss2.subs{end})));
            tmpdata = subsref(tmpMMO.Data.x, ss3);
            fwrite(fid, tmpdata, 'float');
        end
        fclose(fid);
        fprintf('Warning: memory mapped object writing might not be up to date in cache on network drive');
    end
    
    % delete file if necessary
    if ncopies == 1 && obj.writable
        delete(obj.dataFile);
    end
    
    if obj.debug, disp('new copy created, old deleted (length 1)'); end
    obj.dimensions = newdim2;
    obj.dataFile = newFileName;
    obj.writable = true;
    obj = updateWorkspace(obj);
    clear tmpMMO;
    return;
    
else
    % check size to see if it increases
    % ---------------------------------
    newdim1 = obj.dimensions;
    newdim2 = newdim1;
    if length(ss(1).subs) == 1
        if ~ischar(ss(1).subs{1}) && max(ss(1).subs{1}) > prod(newdim1)
            if length(newdim1) > 2
                error('Attempt to grow array along ambiguous dimension.');
            end
        end
        if max(ss(1).subs{1}) > prod(newdim2)
            notOneDim = find(newdim2 > 1);
            if length(notOneDim) == 1
                newdim2(notOneDim) = max(ss(1).subs{1});
            end
        end
    else
        if length(newdim1) == 3 && newdim1(3) == 1, newdim1(end) = []; end
        if length(ss(1).subs) == 2 && length(newdim1) == 3
            if ~ischar(ss(1).subs{2}) && max(ss(1).subs{2}) > prod(newdim1(2:end))
                error('Attempt to grow array along ambiguous dimension.');
            end
            if isnumeric(ss(1).subs{1}), newdim2(1) = max(max(ss(1).subs{1}), newdim2(1)); end
        else
            for index = 1:length(ss(1).subs)
                if isnumeric(ss(1).subs{index})
                    if index > length(newdim2)
                        newdim2(index) = max(ss(1).subs{index});
                    else newdim2(index) = max(max(ss(1).subs{index}), newdim2(index));
                    end
                end
            end
        end
    end
    
    % create new array of new size if necessary
    % -----------------------------------------
    if ~isequal(newdim1, newdim2)
        tmpMMO = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });
        newFileName = mmo.getnewfilename;
        
        % copy file row by row
        % --------------------
        fid = fopen(newFileName, 'w');
        dim1 = [ newdim1 1 1 1 1 ];
        dim2 = [ newdim2 1 1 1 1 ];
        tmpdata1 = zeros(prod(dim2(1:1)) - prod(dim1(1:1)), 1, 'single');
        tmpdata2 = zeros((dim2(2) - dim1(2))*dim2(1), 1, 'single');
        tmpdata3 = zeros((dim2(3) - dim1(3))*prod(dim2(1:2)), 1, 'single');
        
        % copy new data (copy first array)
        % -------------
        for index3 = 1:dim1(3)
            if dim1(1) == dim2(1) && dim1(2) == dim2(2)
                fwrite(fid, tmpMMO.Data.x(:,:,index3), 'float');
            else
                for index2 = 1:dim1(2)
                    if dim1(1) == dim2(1)
                        fwrite(fid, tmpMMO.Data.x(:,index2,index3), 'float');
                    else
                        for index1 = 1:dim1(1)
                            fwrite(fid, tmpMMO.Data.x(index1,index2,index3), 'float');
                        end
                    end
                    fwrite(fid, tmpdata1, 'float');
                end
            end
            fwrite(fid, tmpdata2, 'float');
        end
        fwrite(fid, tmpdata3, 'float');
        fclose(fid);
        fprintf('Warning: memory mapped object writing might not be up to date in cache on network drive');
        
        % delete file if necessary
        if ncopies == 1 && obj.writable
            delete(obj.dataFile);
        end
        
        if obj.debug, disp('new copy created, old deleted'); end
        obj.dimensions = newdim2;
        obj.dataFile = newFileName;
        obj.writable = true;
        clear tmpMMO;
        
        % create copy if necessary
        % ------------------------
    elseif ncopies > 1 || ~obj.writable
        newFileName = mmo.getnewfilename;
        copyfile(obj.dataFile, newFileName);
        obj.dataFile = newFileName;
        obj.writable = true;
        if obj.debug, disp('new copy created'); end
    else
        if obj.debug, disp('using same copy'); end
    end
    
    % copy new data
    tmpMMO = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });
    if ~isa(val, 'mmo')
        tmpMMO.Data.x = builtin('subsasgn', tmpMMO.Data.x, ss, val);
    else
        % copy memory mapped array
        if ndims(val) == 2 && (size(val,1) == 1 || size(val,2) == 1)
            % vector direct copy
            ss2.type = '()';
            ss2.subs = { ':' ':' ':' };
            tmpMMO.Data.x = builtin('subsasgn', tmpMMO.Data.x, ss, subsref(val,ss2));
        else
            ss2.type = '()';
            ss2.subs = { ':' ':' ':' };
            ss3 = ss;
            % array, copy each channel
            for index1 = 1:size(val,1)
                ss2(1).subs{1} = index1;
                if ischar(ss(1).subs{1}) ss3(1).subs{1} = index1;
                else                    ss3(1).subs{1} = ss(1).subs{1}(index1);
                end
                tmpMMO.Data.x = builtin('subsasgn', tmpMMO.Data.x, ss3, subsref(val,ss2));
            end
        end
    end
        
    obj = updateWorkspace(obj);
    
end

