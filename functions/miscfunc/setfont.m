% setfont() - Change all the fonts properties of a figure.
%
% Usage:
%   >>  newdata = setfont( handle, 'key', 'val');
%   >>  [newdata chlab] = setfont( handle, 'key' , 'val', ... );
%   >>  [newdata chlab] = setfont( handle, 'handletype', handletypevalue, 'key' , 'val', ... );
%
% Inputs:
%   handle       - [gcf,gca] figure or plot handle
%   'handletype' - ['xlabels'|'ylabels'|'titles'|'axis'|'strings'] only apply
%                formatting to selected category. Note that this has to be the 
%                first optional input.
%   properties - 'key', 'val' properties for the figure
%
% Example:
%  setfont(gcf, 'fontweight', 'normal', 'fontsize', 14);
%
% Author: Arnaud Delorme, CNL / Salk Institute - SCCN, 25 Oct 2002

%  Copyright (C) 2003 Arnaud Delorme, CNL / Salk Institute - SCCN, La Jolla
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

function setfont(fig, varargin);
    
    if nargin < 1
        help setfont;
        return;
    end
    
    if strcmpi(varargin{1}, 'handletype')
        label = varargin{2};
        varargin = varargin(3:end);
    else
        label = '';
    end
    [hx, hy, hti, hgca, hstr] = findallobjects(fig);
    
    % select a specified category
    % ---------------------------
    if isempty(label)
        h = [hx, hy, hti, hgca, hstr];
    else
        switch lower(label),
         case 'xlabels', h =hx;
         case 'ylabels', h =hy;
         case 'titles',  h =hti;
         case 'axis',    h =hgca;
         case 'strings', h =hstr;
         otherwise, error('Unrecognized ''labels'''); 
        end
    end
    
    % apply formatting
    % ---------------
    for index = 1:length(h)
        isaxis = 0;
        try, get(h(index), 'xtick');  isaxis = 1; catch, end
        if isaxis 
            set(h(index), 'XTickLabelMode', 'manual', 'XTickMode', 'manual');
            set(h(index), 'YTickLabelMode', 'manual', 'YTickMode', 'manual');
        end
        for tmpprop = 1:2:length(varargin)
            if strcmpi(varargin{tmpprop}, 'color') && isaxis
                set(h(index), 'xcolor', varargin{tmpprop+1}, ...
                              'ycolor', varargin{tmpprop+1}, ...
                              'zcolor', varargin{tmpprop+1});                
            else
                try, 
                    set(h(index), varargin{tmpprop}, varargin{tmpprop+1});
                catch, end
            end
        end
    end
    
function [hx, hy, hti, hgca, hstr] = findallobjects(fig);
    handles = findobj(fig)';
    hx   = [];
    hy   = [];
    hti  = [];
    hgca = [];
    hstr = [];
    for index = 1:length(handles)
        try, hx   = [ hx    get(handles(index), 'xlabel')  ]; catch, end
        try, hy   = [ hy    get(handles(index), 'ylabel')  ]; catch, end
        try, hti  = [ hti   get(handles(index), 'title')   ]; catch, end
        try, get(handles(index), 'xtick');  hgca = [ hgca  handles(index) ]; catch, end
        try, get(handles(index), 'string'); hstr = [ hstr  handles(index) ]; catch, end
    end;    
