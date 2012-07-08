% eeglab2loreta() - export chanlocs EEGLAB channel structure
%                   to Tailarach coordinates that can be read using
%                   the LORETA software.
% Usage:
%   eeglab2loreta( chanlocs, data, 'key', 'val', ... );
%
% Inputs:
%   chanlocs       - EEGLAB channel location structure using MNI
%                    coordinates (if not use parameter 'transform').
%   data           - [chan x time] or ICA inverse matrix (EEG.icawinv)
%                    containing component to export. Enter an empty array
%                    to export only the channel location file. Note that
%                    exporting [ chan x time ] requires to set the 'exporterp'
%                    flag.
%
% Optional parameter:
%   'fileloc'      - [string] file name for channel locations. Extension 
%                    '.xyz' (exporting coordinates) or '.txt' (exporting
%                    labels only) will be added automatically. Default is
%                    'loreta_chanlocs'.
%   'filecomp'     - [string] base file name for components. Default is
%                    'comp' so files are written as 'compX.txt' (X representing 
%                    the component number) or 'erp.txt' when exporting ERP.
%   'exporterp'    - ['on'|'off'] export ERP instead of 
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
% Output files:
%   location_file   - contain channel labels or 3-D coordinates in Talairach
%                     space that can be directly read by LORETA
%   component_files - contains 10 rows containing 10 identical component
%                     scalp maps.
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2005
%
% See also: eeglab()

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

function eeglab2loreta( chanlocs, winv, varargin );

    if nargin < 2
        help eeglab2loreta;
        return;
    end;
    
    g = finputcheck( varargin, { 'fileloc'   'string'     []      'loreta_chanlocs';
                                 'filecomp'  'string'     []      'comp';
                                 'compnum'   'integer'    [1 Inf] [];
                                 'transform' 'real'       []      [];
                                 'labelonly' 'string'     { 'on' 'off' } 'off';
                                 'exporterp' 'string'     { 'on' 'off' } 'off';
                                 'excludechan' 'integer'  [1 Inf] [];
                                 'lowchanlim'  'float'    []      NaN });
    if isstr(g), error(g); end;
    
    % remove channels
    % ---------------
    inds = find(cellfun('isempty', { chanlocs.X }));
    g.excludechan = union(g.excludechan, inds);
    if ~isempty(g.excludechan)
        chanlocs(g.excludechan) = [];
        winv(g.excludechan,:)  = [];
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
            chanlocs(rmelec) = [];
            XYZ(rmelec,:)    = [];
            winv(rmelec,:)   = [];
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
    if ~isempty(winv)
        if strcmpi(g.exporterp, 'off')
            if isempty(g.compnum), g.compnum = [1:size(winv,2)]; end;

            for i=g.compnum(:)'
                toLORETA = winv(:,i)'; 
                toLORETA(2:10,:)  = ones(9,1)*toLORETA(1,:);
                filename = sprintf( [ g.filecomp '%d.txt' ], i);
                save('-ascii', filename, 'toLORETA');
            end;
        else
            if strcmpi(g.filecomp, 'comp'), g.filecomp = 'erp'; end;
            toLORETA = double(winv)';
            filename = g.filecomp;
            if ~strcmpi(g.filecomp(end-2:end), 'txt'), filename = [ filename '.txt' ]; end;
            save('-ascii', filename, 'toLORETA');
        end;
    end;
