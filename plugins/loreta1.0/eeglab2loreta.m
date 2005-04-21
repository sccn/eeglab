% eeglab2loreta() - export chanlocs EEGLAB channel structure
%                        to Tailarach coordinates that can be read using
%                        the LORETA software.
% Usage:
%   eeglab2loreta( chanlocs, comps );
%
% Inputs:
%   chanlocs       - EEGLAB channel location structure using MNI
%                    coordinates (if not use parameter 'transform').
%   comps          - ICA inverse matrix or subject of ICA inverse matrix
%                    containing component to export. Enter an empty array
%                    to export only the channel location file.
%
% Optional parameter:
%   'fileloc'      - [string] file name for channel locations. Extension 
%                    '.xyz' (exporting coordinates) or '.txt' (exporting
%                    labels only) will be added automatically. Default is
%                    'loreta_chanlocs'.
%   'filecomp'     - [string] base file name for components. Default is
%                    'compX.txt' (X representing the component number).
%   'compnum'      - [integer array] only export these components.
%   'labelonly'    - ['on'|'off'] only export channel labels and have LORETA
%                    lookup talairach position for them. Default: 'off'.
%   'transform'    - [4x4 matrix] optional homogenous transfromation matrix
%                    to convert channel electrode location to MNI
%                    coordinates.
%   'excludechan'  - [integer] list of channel to omit.
%   'lowchanlim'   - [float] lower z limit (in MNI space) for exporting
%                    channel. For instance 0 will only export channel above
%                    the midline.
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2005
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Arnaud Delorme, SCCN, INC, UCSD, 2005 arno@salk.edu
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

function eeglab2loreta( chanlocs, comps, varargin );

    if nargin < 2
        help chanlocs2talairach;
        return;
    end;
    
    g = finputcheck( varargin, { 'fileloc'   'string'     []      'loreta_chanlocs';
                                 'filecomp'  'string'     []      'comp';
                                 'compnum'   'integer'    [1 Inf] [];
                                 'transform' 'real'       []      [];
                                 'labelonly' 'string'     { 'on' 'off' } [];
                                 'excludechan' 'integer'  [1 Inf] [];
                                 'lowchanlim'  'float'    []      NaN });
    if isstr(g), error(g); end;
    
    % remove channels
    % ---------------
    if ~isemtpy(g.excludechan)
        chanlocs(g.excludechan) = [];
        comps(g.excludechan,:)  = [];
    end;
    
    if strcmpi(g.labelonly, 'off')

        % convert to MNI coordinates
        % --------------------------
        XYZ = [ [chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]'  ];
        if ~isempty(g.transform)
            XYZ = g.transform*[ XYZ ones(size(XYZ,1),1) ]';
            XYZ = XYZ(1:3,:)';
        end;

        % remove all channels below limit
        % -------------------------------
        if ~isnan(g.lowchanlim)
            rmelec = find(XYZ(:,3) < g.lowchanlim);
            XYZ(rmelec,:)   = [];
            comps(rmelec,:) = [];
        end;
        
        % convert to tailairach coordinates
        % ---------------------------------
        talXYZ = mni2tal(XYZ);

        % writing file
        % ------------
        fid = fopen( [ g.fileloc '.xyz' ], 'w');
        fprintf(fid, '%d\n', length(chanlocs));
        for index = 1:length(chanlocs)
            fprintf(fid, '%e %e %e %s\n', talXYZ(index,1), talXYZ(index,2), ...
                    talXYZ(index,3),chanlocs(index).labels);
        end;
        fclose(fid);
        
    else
        
        % export labels only
        % ------------------
        fid = fopen( [ g.fileloc '.txt' ], 'w');
        for index = 1:length(chanlocs)
            fprintf(fid, '%s\n', chanlocs(index).labels);
        end;
        fclose(fid);
        
    end;

    % export components
    % -----------------
    if ~isempty(comps)
        if isempty(g.compnum), g.compnum = [1:size(winv,1)]; end;
            
        for i=g.compnum(:)'
            toLORETA = winv(:,i); 
            filename = sprintf( [ g.filecomp '%d.txt' ], i);
            save('-ascii', filename, 'toLORETA');
        end;
    end;
