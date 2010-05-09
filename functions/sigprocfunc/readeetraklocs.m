% readeetraklocs() - read 3-D location files saved using the EETrak
%                    digitizing software.
% Usage:
%   >> CHANLOCS = readeetraklocs( filename );
%
% Inputs:
%   filename       - [string] file name
%
% Outputs:
%   CHANLOCS       - EEGLAB channel location data structure. 
%                    See help readlocs()
%
% Author: Arnaud Delorme, CNL / Salk Institute, Nov 2003
%
% See also: readlocs()

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function chanlocs = readeetraklocs( filename )
    
    if nargin < 1
        help readeetraklocs;
        return;
    end;
    
    % read location file
    % ------------------
    locs  = loadtxt( filename );
        
    % get label names
    % ---------------
    indlabels = [];
    indpos    = [];
    for ind = 1:size(locs,1)
        if isstr(locs{ind,1}) 
            if strcmpi(locs{ind,1}, 'Labels')
                indlabels = ind;
            end;
            if strcmpi(locs{ind,1}, 'Positions')
                indpos = ind;
            end;
        end;
    end;
    if isempty(indpos) | isempty(indlabels)
        error('Could not find ''Labels'' or ''Position'' tag in electrode file');
    end;
    
    % get positions
    % -------------
    positions = locs(indpos+1:indlabels-1,1:3);
    labels    = locs(indlabels+1:end,:);
        
    % create structure
    % ----------------
    for index = 1:length(labels)
        chanlocs(index).labels = labels{index};
        chanlocs(index).X      = positions{index,1};
        chanlocs(index).Y      = positions{index,2};
        chanlocs(index).Z      = positions{index,3};
    end;
        
    chanlocs = convertlocs(chanlocs, 'cart2all');
