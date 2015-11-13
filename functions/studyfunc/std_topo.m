% std_topo() - uses topoplot() to get the interpolated Cartesian grid of the 
%               specified component topo maps. The topo map grids are saved
%               into a (.icatopo) file and a pointer to the file is stored 
%               in the EEG structure. If such a file already exists, 
%               loads the information from it. 
%
%               Returns the topo map grids of all the requested components. Also
%               returns the EEG sub-structure etc (i.e EEG.etc), which is modified 
%               with a pointer to the float file and some information about the file. 
% Usage:
%               >> X = std_topo(EEG, components, option);  
%
%                  % Returns the ICA topo map grid for a dataset. 
%                  % Updates the EEG structure in the Matlab environment and re-saves
% Inputs:
%   EEG        - an EEG dataset structure. 
%   components - [numeric vector] components in the EEG structure to compute topo maps
%                      {default|[] -> all}      
%   option     - ['gradient'|'laplacian'|'none'] compute gradient or laplacian of
%                the scale topography. This does not acffect the saved file which is
%                always 'none' {default is 'none' = the interpolated topo map}
% Optional inputs
%   'recompute'  - ['on'|'off'] force recomputing topo file even if it is 
%                  already on disk.
%   'fileout'    - [string] Path of the folder to save output. The default
%                  is EEG.filepath
% Outputs:
%   X          - the topo map grid of the requested ICA components, each grid is 
%                     one ROW of X. 
%
% File output: [dataset_name].icatopo
%  
% Authors:  Hilit Serby, Arnaud Delorme, SCCN, INC, UCSD, January, 2005
%
% See also  topoplot(), std_erp(), std_ersp(), std_spec(), std_preclust()

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

function [X] = std_topo(EEG, comps, option, varargin)

if nargin < 1
    help std_topo;
    return;
end;
if isfield(EEG,'icaweights')
   numc = size(EEG.icaweights,1);
else
   error('EEG.icaweights not found');
end
if nargin < 2
   comps = 1:numc;
elseif isempty(comps)
   comps = 1:numc;
end

if nargin < 3
    option = 'none';
end;

g = finputcheck( varargin, { 'recompute'   'string'   { 'on','off' }   'off' ;...
                             'fileout'     'string'   []                EEG.filepath},...
                             'std_topo');
% if isstr(g), error(g); end;

% figure; toporeplot(grid,'style', 'both','plotrad', 0.5, 'intrad', 0.5, 'xsurface' ,Xi, 'ysurface',Yi );

% Topo information found in dataset
% ---------------------------------
if exist(fullfile(g.fileout, [ EEG.filename(1:end-3) 'icatopo' ])) && strcmpi(g.recompute, 'off')
    for k = 1:length(comps)
        tmp = std_readtopo( EEG, 1, comps(k));
        if strcmpi(option, 'gradient')
            [tmpx, tmpy]  = gradient(tmp); %Gradient
            tmp = [tmpx(:); tmpy(:)]';
        elseif strcmpi(option, 'laplacian')
            tmp = del2(tmp); %Laplacian
            tmp = tmp(:)';
        else
            tmp = tmp(:)';
        end;
        
        tmp = tmp(find(~isnan(tmp)));
        if k == 1
            X = zeros(length(comps),length(tmp)) ;
        end
        X(k,:) =  tmp;
    end
    return
end
 
all_topos = [];
for k = 1:numc

    % compute topo map grid (topoimage)
    % ---------------------------------
    chanlocs = EEG.chanlocs(EEG.icachansind);
    if isempty( [ chanlocs.theta ] )
        error('Channel locations are required for computing scalp topographies');
    end;
    [hfig grid plotrad Xi Yi] = topoplot( EEG.icawinv(:,k), chanlocs, ...
                                          'verbose', 'off',...
                                           'electrodes', 'on' ,'style','both',...
                                           'plotrad',0.55,'intrad',0.55,...
                                           'noplot', 'on', 'chaninfo', EEG.chaninfo);

    all_topos = setfield(all_topos, [ 'comp' int2str(k) '_grid' ], grid);
    all_topos = setfield(all_topos, [ 'comp' int2str(k) '_x' ]   , Xi(:,1));
    all_topos = setfield(all_topos, [ 'comp' int2str(k) '_y' ]   , Yi(:,1));
    
end

% Save topos in file
% ------------------
all_topos.datatype = 'TOPO';
tmpfile = fullfile( g.fileout, [ EEG.filename(1:end-3) 'icatopo' ]); 
std_savedat(tmpfile, all_topos);

for k = 1:length(comps)
    tmp =  getfield(all_topos, [ 'comp' int2str(comps(k)) '_grid' ]);
    
    if strcmpi(option, 'gradient')
        [tmpx, tmpy]  = gradient(tmp); % Gradient
        tmp = [tmpx(:); tmpy(:)]';
    elseif strcmpi(option, 'laplacian')
        tmp = del2(tmp); % Laplacian
        tmp = tmp(:)';
    else
        tmp = tmp(:)';
    end;

    tmp = tmp(find(~isnan(tmp)));
    if k == 1
        X = zeros(length(comps),length(tmp)) ;
    end
    X(k,:) =  tmp;
end
