% readneurolocs() - read neuroscan polar location file (.asc)
%
% Usage:
%   >> CHANLOCS = readneurolocs( filename );
%   >> CHANLOCS = readneurolocs( filename, 'key1', val1, ...);
%
% Inputs:
%   filename       - file name or matlab cell array { names x_coord y_coord }
%
% Optional inputs:
%   same as caliblocs()
%   note that if no optional input are provided, re-centering will be
%   performed automatically and re-scaling of coordinates will be
%   performed for '.asc' files (not '.dat' files).
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
% Revision 1.10  2003/12/02 19:15:48  arno
% pass on paramters to adjustlocs
%
% Revision 1.9  2003/12/02 02:34:10  arno
% now using adjustlocs
%
% Revision 1.8  2003/12/01 17:10:17  arno
% automatic calibration
%
% Revision 1.7  2003/12/01 02:31:01  arno
% correct slight innacuracy when reading file
%
% Revision 1.6  2003/08/05 18:24:04  arno
% header
%
% Revision 1.5  2003/07/29 21:25:20  arno
% debug C3-C4
%
% Revision 1.4  2003/07/29 21:05:53  arno
% handle cell array for channel locations
%
% Revision 1.3  2003/06/13 00:53:47  arno
% debuging
%
% Revision 1.2  2003/03/04 20:04:17  arno
% debuging
%
% Revision 1.1  2003/03/04 19:18:35  arno
% Initial revision
%

function chanlocs = readneurolocs( filename, varargin)
    
    if nargin < 1
        help readneurolocs;
        return;
    end;
    if nargin < 2
        plottag = 0;
    end;
    
    % read location file
    % ------------------
    if isstr( filename )
        locs  = loadtxt( filename );
    
        if locs{1,1}(1) == ';'
            % remove trailing control channels
            % --------------------------------
            while isnumeric( locs{end,1} ) & locs{end,1} ~= 0
                locs  = locs(1:end-1,:);
            end;    
            
            % find first numerical index
            % --------------------------
            index = 1;
            while isstr( locs{index,1} )
                index = index + 1;
            end;
            
            % extract location array
            % ----------------------   
            nchans = size( locs, 1 ) - index +1;
            chans  = cell2mat(locs(end-nchans+1:end, 1:5));
            names  = locs(end-nchans*2+1: end-nchans, 2);
            for index = 1:length(names)
                if ~isstr(names{index})
                    names{index} = int2str(names{index});
                end;
            end;
            x = chans(:,3);
            y = -chans(:,4);
        else
            [tmp2 tmpind] = sort(cell2mat(locs(:,1))');
            locs = locs(tmpind,:);
            y      = cell2mat(locs(:,end));
            x      = cell2mat(locs(:,end-1));
            x      = x/513.1617*44;
            y      = y/513.1617*44;
            names = locs(:,end-2);
        end;
    else
        names = filename{1};
        x     = filename{2};
        y     = filename{3};
    end;
    
    % second solution using angle
    % ---------------------------
    [phi,theta] = cart2pol(x, y);
    phi = phi/pi*180;
    
    % convert to other types of coordinates
    % -------------------------------------
    labels = names';
    chanlocs = struct('labels', labels, 'sph_theta_besa', mat2cell(theta)', 'sph_phi_besa', mat2cell(phi)');
    chanlocs = convertlocs( chanlocs, 'sphbesa2all');
    for index = 1:length(chanlocs)
        chanlocs(index).labels = num2str(chanlocs(index).labels);
    end;
    
    % re-calibration
    % --------------
    chanlocs = adjustlocs(chanlocs, 'autoscale', 'on', 'autorotate', 'off', varargin{:});
