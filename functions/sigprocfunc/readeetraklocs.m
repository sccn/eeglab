% readeetraklocs() - read 3-d location files saved using the EETrak
%                    digitizing software.
%
% Usage:
%   >> CHANLOCS = readeetraklocs( filename );
%
% Inputs:
%   filename       - [string] file name
%
% Outputs:
%   CHANLOCS       - EEGLAB channel location data structure. See
%                    help readlocs()
%
% Author: Arnaud Delorme, CNL / Salk Institute, Nov 2003
%
% See also: readlocs()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $

function chanlocs = readeetraklocs( filename)
    
    if nargin < 1
        help readeetraklocs;
        return;
    end;
    
    % read location file
    % ------------------
    locs  = loadtxt( filename );
        
    % get label names
    % ---------------
    labels = locs(end,:);
    
    % get positions
    % -------------
    positions = locs(4:end-2,1:3);
        
    % create structure
    % ----------------
    for index = 1:length(labels)
        chanlocs(index).labels = labels{index};
        chanlocs(index).X      = positions{index,1};
        chanlocs(index).Y      = positions{index,2};
        chanlocs(index).Z      = positions{index,3};
    end;
        
    chanlocs = convertlocs(chanlocs, 'cart2all');