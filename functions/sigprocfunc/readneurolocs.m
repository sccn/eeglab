% readneurolocs() - read neuroscan polar location file (.asc)
%
% Usage:
%   >> CHANLOCS = readneurolocs( filename );
%   >> CHANLOCS = readneurolocs( filename, plot, czindex, fzindex );
%
% Inputs:
%   filename       - file name
%
% Optional inputs:
%   plot           - [0|1] if 1, plot the electrode locations
%   czindex        - index of electrode Cz if the label is not defined 
%                    in the file
%   fzindex        - index of electrode Fz if the label is not defined 
%                    in the file. (for 3D accurate conversion of electrode 
%                    positions)
%
% Outputs:
%   CHANLOCS       - EEGLAB channel location data structure. See
%                    help readlocs()
%
% Author: Arnaud Delorme, CNL / Salk Institute, 4 March 2003
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

function readneurolocs( filename, plottag, indexcz, indexfz)
    
    if nargin < 1
        help readneurolocs;
        return;
    end;
    if nargin < 2
        plottag = 0;
    end;
    
    % read location file
    % ------------------
    locs = loadtxt( filename );
    
    % find first numerical index
    % --------------------------
    index = 1;
    while isstr( locs{index,1} )
        index = index + 1;
    end;
        
    % extract location array
    % ----------------------   
    nchans = size( locs, 1 ) - index +1;
    chans  = locs(end-nchans+1:end, 1:5);
    names  = locs(end-nchans*2+1: end-nchans, 2);
    
    % find Cz and Fz
    % --------------
    if exist('czindex') == 1
        indexcz = strmatch( 'cz', lower(names') );
    end;
    if exist('czindex') == 1
        indexfz = strmatch( 'fz', lower(names') );
    end;
    
    % plot all channels
    % -----------------
    x = chans(:,3);
    y = chans(:,4);    
    if plottag
        figure;
        for index = 1:length(x)
            plot( x(index), y(index), '+');
            hold on; % do not erase
            text( x(index)+0.01, y(index), int2str(index));
        end;
    end;
    
    % convert to polar
    % ----------------
    x = x - x(indexcz);
    y = - y + y(indexfz);    
    [theta r] = cart2pol(x, y);
    
    % skink radius 
    % ------------
    r = r/r(indexfz)*0.25556
    
    % create structure
    % ----------------
    for index = 1:length(names)
        chanlocs(index).labels = names{index};
        chanlocs(index).theta  = theta(index)/pi*180+90;
        chanlocs(index).radius = r(index);
    end;
        