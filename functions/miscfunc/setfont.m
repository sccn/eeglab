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
%                formating to selected category. Note that this has to be the 
%                first optional input.
%   properties - 'key', 'val' properties for the figure
%
% Exemple:
%  setfont(gcf, 'fontweight', 'normal', 'fontsize', 14);
%
% Author: Arnaud Delorme, CNL / Salk Institute - SCCN, 25 Oct 2002

%  Copyright (C) 2003 Arnaud Delorme, CNL / Salk Institute - SCCN, La Jolla
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

function setfont(fig, varargin);
    
    if nargin < 1
        help setfont;
        return;
    end;
    
    if strcmpi(varargin{1}, 'handletype')
        label = varargin{2};
        varargin = varargin(3:end);
    else
        label = '';
    end;
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
        end;
    end;
    
    % apply formating
    % ---------------
    for index = 1:length(h)
        isaxis = 0;
        try, get(h(index), 'xtick');  isaxis = 1; catch, end;
        if isaxis 
            set(h(index), 'XTickLabelMode', 'manual', 'XTickMode', 'manual');
            set(h(index), 'YTickLabelMode', 'manual', 'YTickMode', 'manual');
        end;
        for tmpprop = 1:2:length(varargin)
            if strcmpi(varargin{tmpprop}, 'color') & isaxis
                set(h(index), 'xcolor', varargin{tmpprop+1}, ...
                              'ycolor', varargin{tmpprop+1}, ...
                              'zcolor', varargin{tmpprop+1});                
            else
                try, 
                    set(h(index), varargin{tmpprop}, varargin{tmpprop+1});
                catch, end;
            end;
        end;
    end;
    
function [hx, hy, hti, hgca, hstr] = findallobjects(fig);
    handles = findobj(fig)';
    hx   = [];
    hy   = [];
    hti  = [];
    hgca = [];
    hstr = [];
    for index = 1:length(handles)
        try, hx   = [ hx    get(handles(index), 'xlabel')  ]; catch, end;
        try, hy   = [ hy    get(handles(index), 'ylabel')  ]; catch, end;
        try, hti  = [ hti   get(handles(index), 'title')   ]; catch, end;
        try, get(handles(index), 'xtick');  hgca = [ hgca  handles(index) ]; catch, end;
        try, get(handles(index), 'string'); hstr = [ hstr  handles(index) ]; catch, end;
    end;    
