% EEGLAB_EXECMENU - function called to execute EEGLAB menu
%
% Usage:
%   eeglab_execmenu(label, func, parameter);
%
% Inputs:
%   label     - [string] label of the menu
%   func      - [string] name of the pop function. The function extract the
%               callback from the menu, find the function and add the
%               parameters below before executing the callback.
%   parameter - [cell] parameter to the function (in addition to existing
%               parameters.
%
% Output: the command is executed in the global workspace. Usually the 
%         EEG structure is affected.
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2022
%
% See also: EEGLAB

% Copyright (C) 2022 Arnaud Delorme
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
 
function eeglab_execmenu(label, func, parameters)

if nargin < 3
    fprintf(2, 'This function requires 3 inputs\n')
    help eeglab_execmenu;
    return
end

% find figure
W_MAIN = findobj('tag', 'EEGLAB');
if isempty(W_MAIN)
    error('EEGLAB figure not found');
end

% find menu
menuhandle = findobj(W_MAIN, 'Label', label);
if isempty(menuhandle)
    error('Could not find menu with label ''%s''', label);
end

% check if menu contains command
cb = get(menuhandle, 'callback');
if iscell(cb) 
    if contains(cb{1}, func)
        cb = cb{1};
    else
        cb = cb{end};
    end
end
ind = strfind(cb, func);
if length(ind) ~= 1
    fprintf(2, 'Menu with label ''%s'' found, but the callback does not contain the string ''%s''\n', label, func);
    ind = ind(end);
end

% extract callback and add parameters
if ~isempty(parameters)
    if cb(ind+length(func)) == ';'
        cb = [ cb(1:ind+length(func)-1) '()' cb(ind+length(func):end) ];
    end
    indEnd = find(cb(ind:end) == ')'); 
    indEnd = indEnd(1);
    if cb(ind+indEnd-2) == '(' % no parameter before
        addDelim = '';
    else
        addDelim = ',';
    end
    cb = [cb(1:ind+indEnd-2) addDelim vararg2str(parameters) cb(ind+indEnd-1) cb(ind+indEnd:end)];
end

% execute in global workspace (as menu would)
evalin('base', 'DEBUG_EEGLAB_MENUS = true;');
evalin('base', cb);
evalin('base', 'clear DEBUG_EEGLAB_MENUS;');
