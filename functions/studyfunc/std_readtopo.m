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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [X, yi, xi ] = std_readtopo(ALLEEG, abset, comps, option, mode)

X = [];
yi = [];
xi = [];
if nargin < 4
    option = 'none';
end
if nargin < 5
    mode = '2Dmap';
end
filename = correctfile(fullfile( ALLEEG(abset).filepath,[ ALLEEG(abset).filename(1:end-3) 'icatopo']),ALLEEG(abset).filepath);
tmpfile  = which(filename);
if ~isempty(tmpfile), filename = tmpfile; end

% 061411, 2:51pm
% Modified by Joaquin
% while getfield(dir(filename), 'bytes') < 1000
i = 1;
while getfield(dir(filename), 'bytes') < 5000
    topo = load( '-mat', filename);
    filename = correctfile(topo.file, ALLEEG(abset).filepath);
    tmpfile  = which(filename);
    if ~isempty(tmpfile), filename = tmpfile; end
    if(i>100) 
        error('too many attempts to find valid icatopo');
    end
    i = i+1;
end

for k = 1:length(comps)

    if length(comps) < 3
        try
            topo = load( '-mat', filename, ...
                         [ 'comp' int2str(comps(k)) '_grid'], ...
                         [ 'comp' int2str(comps(k)) '_x'], ...
                         [ 'comp' int2str(comps(k)) '_y'] );
        catch
            error( [ 'Cannot read file ''' filename '''' ]);
        end
    elseif k == 1
        try
            topo = load( '-mat', filename);
        catch
            error([ 'Missing scalp topography file - also necessary for ERP polarity' 10 'Try recomputing scalp topographies for components' ]);
        end
    end
    
    try,
        tmp =  getfield(topo, [ 'comp' int2str(comps(k)) '_grid' ]);
    catch,
        error([ 'Empty scalp topography file - also necessary for ERP polarity' 10 'Try recomputing scalp topographies for components' ]);
    end
        
    if strcmpi(option, 'gradient')
        [tmpx, tmpy]  = gradient(tmp); % Gradient
        tmp        = tmpx;
        tmp(:,:,2) = tmpy;
    elseif strcmpi(option, 'laplacian')
        tmp = del2(tmp); % Laplacian
    end

    if length(comps) > 1 || strcmpi(mode, 'preclust')
        tmp = tmp(find(~isnan(tmp))); % remove NaN for more than 1 component
    end
    if k == 1
        X = zeros([ length(comps) size(tmp) ]) ;
    end
    X(k,:,:,:) =  tmp;
    if k == 1 
        yi   = getfield(topo, [ 'comp' int2str(comps(k)) '_y']);
        xi   = getfield(topo, [ 'comp' int2str(comps(k)) '_x']);
    end
end
X = squeeze(X);

return;

function filename = correctfile(filename, datasetpath)
    comp = computer;
    if filename(2) == ':' && ~strcmpi(comp(1:2), 'PC') 
        filename = [filesep filename(4:end) ];
        filename(find(filename == '\')) = filesep;
    end
    
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
                    error([ 'Cannot load file ''' [ tmpf ext ] '''' 10 'Go back and recompute the data file.' 10 'Note that plotting ICA component ERPs require' 10 'to precompute ICA topographies (see tutorial)']);
                end
            end
        end
    end;        
    
