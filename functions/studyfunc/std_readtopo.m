% std_readtopo() - returns the scalp map of a specified ICA component, assumed
%                  to have been saved in a Matlab file, [dataset_name].icatopo, 
%                  in the same directory as the dataset file. If this file does 
%                  not exist, use std_topo() to create it, else a pre-clustering 
%                  function that calls it: pop_preclust() or eeg_preclust().  
% Usage:    
%   >> [grid, y, x ] = std_readtopo(ALLEEG, setindx, component);  
%   >> [grid, y, x ] = std_readtopo(ALLEEG, setindx, component, transform, mode);  
%
% Inputs:
%   ALLEEG     - vector of EEG datasets (can also be one EEG set). 
%                must contain the dataset of interest (see 'setindx' below).
%   setindx    - [integer] an index of an EEG dataset in the ALLEEG
%                structure, for which to get the component ERP.
%   component  - [integer] index of the component for which the scalp map 
%                grid should be returned. 
%   transform  - ['none'!'laplacian'|'gradient'] transform scalp map to
%                laplacian or gradient map. Default is 'none'.
%   mode       - ['2dmap'|'preclust'] return either a 2-D array for direct
%                plotting ('2dmap') or an array formated for preclustering
%                with all the NaN values removed (ncomps x points). Default
%                is '2dmap' for 1 component and 'preclust' for several.
%
% Outputs:
%   grid      - square scalp-map color-value grid for the requested ICA component 
%               in the specified dataset, an interpolated Cartesian grid as output 
%               by topoplot(). 
%   y         - y-axis values for the interpolated grid
%   x         - x-axis values of the interpolated grid
%
%  See also  std_topo(), std_preclust()
%
% Authors: Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, February, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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
% Revision 1.14  2007/12/09 01:58:23  arno
% nothing
%
% Revision 1.13  2007/12/09 01:46:59  arno
% allow reading windows files under Linux
%
% Revision 1.12  2007/12/09 01:02:22  arno
% fix reading topo file
%
% Revision 1.11  2007/12/09 00:10:40  arno
% [6~loading the correct file
%
% Revision 1.10  2007/11/14 02:42:48  arno
% fix the file reference thing
%
% Revision 1.9  2007/09/11 10:53:25  arno
% handles now one file per dataset
%
% Revision 1.8  2006/03/14 02:39:40  scott
% help msg
%
% Revision 1.7  2006/03/11 07:23:51  arno
% header
%
% Revision 1.6  2006/03/10 00:37:17  arno
% error msg
%
% Revision 1.5  2006/03/09 19:00:42  arno
% reading Matlab file
%
% Revision 1.4  2006/03/09 00:00:54  arno
%  now saving Matlab file
%
% Revision 1.3  2006/03/08 20:32:48  arno
% rename func
%
% Revision 1.2  2006/03/07 22:14:40  arno
% use fullfile
%

function [X, yi, xi ] = std_readtopo(ALLEEG, abset, comps, option, mode)

X = [];
yi = [];
xi = [];
if nargin < 4
    option = 'none';
end;
if nargin < 5
    mode = '2Dmap';
end;
filename = correctfile(fullfile( ALLEEG(abset).filepath,[ ALLEEG(abset).filename(1:end-3) 'icatopo']));
while (getfield(dir(filename), 'bytes') < 1000)
    topo = load( '-mat', filename);
    filename = correctfile(topo.file);
end;

for k = 1:length(comps)

    if length(comps) < 3
        warning off;
        try
            topo = load( '-mat', filename, ...
                         [ 'comp' int2str(comps(k)) '_grid'], ...
                         [ 'comp' int2str(comps(k)) '_x'], ...
                         [ 'comp' int2str(comps(k)) '_y'] );
        catch
            error( [ 'Cannot read file ''' filename '''' ]);
        end;
        warning on;
    elseif k == 1
        try
            topo = load( '-mat', filename);
        catch
            error( [ 'Cannot read file ''' filename '''' ]);
        end;
    end;
    
    tmp =  getfield(topo, [ 'comp' int2str(comps(k)) '_grid' ]);
    
    if strcmpi(option, 'gradient')
        [tmpx, tmpy]  = gradient(tmp); % Gradient
        tmp        = tmpx;
        tmp(:,:,2) = tmpy;
    elseif strcmpi(option, 'laplacian')
        tmp = del2(tmp); % Laplacian
    end;

    if length(comps) > 1 | strcmpi(mode, 'preclust')
        tmp = tmp(find(~isnan(tmp))); % remove NaN for more than 1 component
    end;
    if k == 1
        X = zeros([ length(comps) size(tmp) ]) ;
    end
    X(k,:,:,:) =  tmp;
    if k == 1 
        yi   = getfield(topo, [ 'comp' int2str(comps(k)) '_y']);
        xi   = getfield(topo, [ 'comp' int2str(comps(k)) '_x']);
    end;
end
X = squeeze(X);

return;

function filename = correctfile(filename)
    if filename(2) == ':' & ~strcmpi(computer, 'pcwin')
        filename = filename(4:end);
        filename(find(filename == '\')) = filesep;
    end;
    if ~exist(filename)
        [tmpp tmpf ext] = fileparts(filename);
        if exist([tmpf ext])
            filename = [tmpf ext];
        else
            [tmpp2 tmpp1] = fileparts(tmpp);
            if exist(fullfile(tmpp1, [ tmpf ext ]))
                filename = fullfile(tmpp1, [ tmpf ext ]);
            else
                error([ 'Cannot load file ''' [ tmpf ext ] '''' ]);
            end;
        end;
    end;        
    
