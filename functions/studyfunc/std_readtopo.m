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
tmpfile  = which(filename);
if ~isempty(tmpfile), filename = tmpfile; end;

% 061411, 2:51pm
% Modified by Joaquin
% while getfield(dir(filename), 'bytes') < 1000
i = 1;
while getfield(dir(filename), 'bytes') < 5000
    topo = load( '-mat', filename);
    filename = correctfile(topo.file, ALLEEG(abset).filepath);
    tmpfile  = which(filename);
    if ~isempty(tmpfile), filename = tmpfile; end;
    if(i>100) 
        error('too many attempts to find valid icatopo');
    end
    i = i+1;
end;

for k = 1:length(comps)

    if length(comps) < 3
        try
            topo = load( '-mat', filename, ...
                         [ 'comp' int2str(comps(k)) '_grid'], ...
                         [ 'comp' int2str(comps(k)) '_x'], ...
                         [ 'comp' int2str(comps(k)) '_y'] );
        catch
            error( [ 'Cannot read file ''' filename '''' ]);
        end;
    elseif k == 1
        try
            topo = load( '-mat', filename);
        catch
            error([ 'Missing scalp topography file - also necessary for ERP polarity' 10 'Try recomputing scalp topographies for components' ]);
        end;
    end;
    
    try,
        tmp =  getfield(topo, [ 'comp' int2str(comps(k)) '_grid' ]);
    catch,
        error([ 'Empty scalp topography file - also necessary for ERP polarity' 10 'Try recomputing scalp topographies for components' ]);
    end;
        
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

function filename = correctfile(filename, datasetpath)
    comp = computer;
    if filename(2) == ':' & ~strcmpi(comp(1:2), 'PC') 
        filename = [filesep filename(4:end) ];
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
                filename = fullfile(datasetpath, [ tmpf ext ]);
                if ~exist(filename)
                    error([ 'Cannot load file ''' [ tmpf ext ] '''' ]);
                end;
            end;
        end;
    end;        
    
