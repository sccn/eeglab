% questdlg2() - questdlg function clone with coloring and help for 
%               eeglab().
%
% Usage: same as questdlg()
%
% Warning: 
% Case of button text and result might be changed by the function
%
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 11 August 2002
%
% See also: inputdlg2(), errordlg2(), supergui(), inputgui()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, arno@salk.edu
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

function [result] = questdlg2(Prompt,Title,varargin);

result = '';
if nargin < 2
   help questdlg2;
   return;
end
if isempty(varargin)
    varargin = { 'Yes' 'No' 'Cancel' 'Yes' };
end
result = varargin{end};
if Prompt(end) == 10, Prompt(end) = []; end
if Prompt(end) == 10, Prompt(end) = []; end
if Prompt(end) == 10, Prompt(end) = []; end
if Prompt(end) == 10, Prompt(end) = []; end

fig = figure('visible', 'off');
set(gcf, 'name', Title);

listui = {};
geometry = {};
if ~isempty(find(Prompt == 10))
    indlines = find(Prompt == 10);
    if indlines(1) ~= 1, indlines = [ 0 indlines ]; end
    if indlines(end) ~= length(Prompt), indlines = [ indlines length(Prompt)+1 ]; end
    for index = 1:length(indlines)-1
        geometry{index} = [1];
        listui{index} = { 'Style', 'text', 'string' Prompt(indlines(index)+1:indlines(index+1)-1) };
    end
else
    for index = 1:size(Prompt,1)
        geometry{index} = [1];
        listui{index} = { 'Style', 'text', 'string' Prompt(index,:) };
    end
end
listui{end+1} = {};

geometry = { geometry{:} 1 ones(1,length(varargin)-1) };
for index = 1:length(varargin)-1 % ignoring default val
	listui = {listui{:} { 'width',80,'align','center','Style', 'pushbutton', 'string', varargin{index}, 'callback', ['set(gcbf, ''userdata'', ''' varargin{index} ''');'] }  };
	if strcmp(varargin{index}, varargin{end})
		listui{end}{end+1} = 'fontweight';
		listui{end}{end+1} = 'bold';
	end
end

%cr = length(find(Prompt == char(10)))+1;
%if cr == 1
%	cr = size(Prompt,1);
%end
%cr = cr^(7/);
%if cr >= 8, cr = cr-1; end
%if cr >= 4, cr = cr-1; end
%[tmp tmp2 allobj] = supergui( 'fig', fig, 'geomhoriz', geometry, 'geomvert', [cr 1 1], 'uilist', listui, ...
[tmp tmp2 allobj] = supergui( 'fig', fig, 'geomhoriz', geometry, 'uilist', listui, ...
    'borders', [0.05 0.015 0.08 0.06], 'spacing', [0 0], 'horizontalalignment', 'left', 'adjustbuttonwidth', 'off' );

waitfor( fig, 'userdata');
try,
	result = get(fig, 'userdata');
	close(fig);
    drawnow;
end
