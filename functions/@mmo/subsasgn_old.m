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

% check stack
% -----------
stack = dbstack;
stack(1) = [];
stack    = rmfield(stack, 'line');
ncopies = 0;
if ~isempty(stack)
    % check if we are in a different workspace
    if ~isequal(stack, obj.workspace)
        % if subfunction, must be a copie
        if ~isempty(obj.workspace) && strcmpi(stack(end).file, obj.workspace(end).file) && ...
                ~strcmpi(stack(end).name, obj.workspace(end).name)
            % We are within a subfunction. The MMO must have
            % been passed as an argument (otherwise the current
            % workspace and the workspace variable would be
            % equal).
            ncopies = 2;
        else
            tmpVar = evalin('caller', 'nargin;'); % this does not work
            if ~isscript(stack(1).file)
                ncopies = 2;
                % we are within a function. The MMO must have
                % been passed as an argument (otherwise the current
                % workspace and the workspace variable would be
                % equal).
            else
                % we cannot be in a function with 0 argument
                % (otherwise the current workspace and the workspace
                % variable would be equal). We must assume that
                % we are in a script.
                while ~isempty(stack) && ~isequal(stack, obj.workspace)
                    stack(1) = [];
                end
                if ~isequal(stack, obj.workspace)
                    ncopies = 2;
                end
            end
        end
    end
end

% check local variables
% ---------------------
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
    newFileName = sprintf('memapdata_%.9d%.9d.fdt', round(rand(1)*10^9), round(rand(1)*10^9));
    
    % find non singleton dimension
    % ----------------------------
    nonSingleton = [];
    for index = 1:length(subs)
        if ~ischar(subs{index}) % can only be ":"
            nonSingleton(end+1) = index;
            subs2 = setdiff_bc([1:newdim1(index)], subs{index}); % invert selection
        end
    end
    if length(nonSingleton) > 1, error('A null assignment can have only one non-colon index'); end
    if isempty(nonSingleton), obj = []; return; end
    
    % compute new final dimension
    % ---------------------------
    if length(ss(1).subs) == 1
        fid = fopen(newFileName, 'w');
        newdim2 = prod(newdim2)-length(ss(1).subs{1});
        newindices = setdiff_bc([1:prod(newdim1)], ss(1).subs{1});
        for index = newindices
            fwrite(fid, tmpMMO.Data.x(index), 'float');
        end
        fclose(fid);
    else
        newdim2(nonSingleton) = newdim2(nonSingleton)-length(subs{index});
        tmpdata = builtin('subsref', tmpMMO.Data.x, s);
        fid = fopen(newFileName, 'w');
        fwrite(fid, tmpMMO.Data.x(index), 'float');
        fclose(fid);
    end
    
    % delete file if necessary
    if ncopies == 1 && obj.writable
        delete(obj.dataFile);
    end
    
    if obj.debug, disp('new copy created, old deleted (length 1)'); end
    obj.dimensions = [1 newdim2];
    obj.dataFile = newFileName;
    obj.writable = true;
    clear tmpMMO;
    return;
    
elseif ncopies == 1 && obj.writable
    for index = 1:length(ss(1).subs)
        newdim2(notOneDim) = newdim2(notOneDim)-length(ss(1).subs{1});
        if index > length(newdim2)
            newdim2(index) = max(ss(1).subs{index});
        else newdim2(index) = max(max(ss(1).subs{index}), newdim2(index));
        end
    end
    
    % create new array of new size
    % -----------------------------------------
    tmpMMO = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });
    newFileName = sprintf('memapdata_%.9d%.9d.fdt', round(rand(1)*10^9), round(rand(1)*10^9));
    
    % copy file row by row
    % --------------------
    fid = fopen(newFileName, 'w');
    tmpdata = zeros(prod(newdim2(1)) - prod(newdim1(1)), 1, 'single');
    for index = 1:prod(newdim1(2:end))
        fwrite(fid, tmpMMO.Data.x(:,index), 'float');
        fwrite(fid, tmpdata, 'float');
    end
    if prod(newadim1(2:end)) ~= prod(newdim2(2:end))
        tmpdata = zeros(newdim2(1), 1, 'single');
        for index = prod(newdim1(2:end))+1:prod(newdim2(2:end))
            fwrite(fid, tmpdata, 'float');
        end
    end
    fclose(fid);
    
    % delete file if necessary
    if ncopies == 1 && obj.writable
        delete(obj.dataFile);
    end
    
    if obj.debug, disp('new copy created, old deleted'); end
    obj.dimensions = newdim2;
    obj.dataFile = newFileName;
    obj.writable = true;
    clear tmpMMO;
    
else
    % check size
    % ----------
    newdim1 = obj.dimensions;
    newdim2 = obj.dimensions;
    if length(ss(1).subs) == 1
        if max(ss(1).subs{1}) > prod(newdim2)
            notOneDim = find(newdim2 > 1);
            if length(notOneDim) == 1
                newdim2(notOneDim) = max(ss(1).subs{1});
            end
        end
    else
        for index = 1:length(ss(1).subs)
            if index > length(newdim2)
                newdim2(index) = max(ss(1).subs{index});
            else newdim2(index) = max(max(ss(1).subs{index}), newdim2(index));
            end
        end
    end
    
    % create new array of new size if necessary
    % -----------------------------------------
    if ~isequal(newdim1, newdim2)
        tmpMMO = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });
        newFileName = sprintf('memapdata_%.9d%.9d.fdt', round(rand(1)*10^9), round(rand(1)*10^9));
        
        % copy file row by row
        % --------------------
        fid = fopen(newFileName, 'w');
        tmpdata = zeros(prod(newdim2(1)) - prod(newdim1(1)), 1, 'single');
        for index = 1:prod(newdim1(2:end))
            fwrite(fid, tmpMMO.Data.x(:,index), 'float');
            fwrite(fid, tmpdata, 'float');
        end
        if prod(newdim1(2:end)) ~= prod(newdim2(2:end))
            tmpdata = zeros(newdim2(1), 1, 'single');
            for index = prod(newdim1(2:end))+1:prod(newdim2(2:end))
                fwrite(fid, tmpdata, 'float');
            end
        end
        fclose(fid);
        
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
        newFileName = sprintf('memapdata_%.9d%.9d.fdt', round(rand(1)*10^9), round(rand(1)*10^9));
        copyfile(obj.dataFile, newFileName);
        obj.dataFile = newFileName;
        obj.writable = true;
        if obj.debug, disp('new copy created'); end
    else
        if obj.debug, disp('using same copy'); end
    end
    
    tmpMMO = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });
    tmpMMO.Data.x = builtin('subsasgn', tmpMMO.Data.x, ss, val);
    
end

return;

%     i.type = '()';
%     i.subs = { ':' ':' ':' };
%     res = subsref(obj, i); % note that subsref cannot be called directly


% subfunction checking the number of copies
% -----------------------------------------
function ncopies = checkcopies_local(obj, arg);
ncopies = 0;
if isstruct(arg)
    for ilen = 1:length(arg)
        for index = fieldnames(arg)'
            ncopies = ncopies + checkcopies_local(obj, arg(ilen).(index{1}));
            if ncopies > 1, return; end
        end
    end
elseif iscell(arg)
    for index = 1:length(arg(:))
        ncopies = ncopies + checkcopies_local(obj, arg{index});
        if ncopies > 1, return; end
    end
elseif isa(arg, 'mmo') && isequal(obj, arg)
    ncopies = 1;
else
    ncopies = 0;
end
