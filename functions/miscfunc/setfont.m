% setfont() - Change all the fonts properties of a figure.
%
% Usage:
%   >>  newdata = setfont( handle, 'key', 'val');
%   >>  [newdata chlab] = setfont( handle, 'key' , 'val', ... );
%   >>  [newdata chlab] = setfont( handle, 'labels', 'key' , 'val', ... );
%
% Inputs:
%   handle     - [gcf,gca] figure or plot handle
%   'labels'   - if the keyword label is present, the formating is only
%                applied to x and y labels and titles.
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

% $Log: not supported by cvs2svn $

function setfont(fig, varargin);
    
    if nargin < 1
        help setfont;
        return;
    end;
    
    if strcmpi(varargin{1}, 'labels')
        varargin = varargin(2:end);
        label = 1;
    else
        label = 0;
    end;
    h = findallobjects(fig, label);
    
    for index = 1:length(h)
        try, 
            set(h(index), 'XTickLabelMode', 'manual', 'XTickMode', 'manual');
            set(h(index), 'YTickLabelMode', 'manual', 'YTickMode', 'manual');
        catch, end;
        try, 
            set(h(index), varargin{:});
        catch, end;
    end;
    
function handles = findallobjects(fig, label);
    handles = findobj(fig)';
    lhandles = [];
    for index = 1:length(handles)
        try, 
            lhandles = [ lhandles get(handles(index), 'xlabel')  ];
            lhandles = [ lhandles get(handles(index), 'ylabel')  ];
            lhandles = [ lhandles get(handles(index), 'title')  ];
        catch, 
        end;
    end;    
    if nargin < 2 | label == 1
        handles = lhandles;
    else 
        handles = [ handles lhandles ];
    end;
