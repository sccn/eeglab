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
%
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
% Authors:  Arnaud Delorme, SCCN, INC, UCSD, January, 2005
%
% See also  topoplot(), std_erp(), std_ersp(), std_spec(), std_preclust()

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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

function [X] = std_topo(EEG, comps, varargin)

if nargin < 1
    help std_topo;
    return;
end
if nargin > 2 || strcmpi(varargin{1}, 'none') || strcmpi(varargin{1}, 'gradient') || strcmpi(varargin{1}, 'laplacian') 
    fprintf('std_topo: option to return gradient or laplacian is deprecated and ignored\n');
    varargin = varargin(2:end);
end

g = finputcheck( varargin, { 'recompute'   'string'   { 'on','off' }   'off' ; ...
                             'trialinfo'   'struct'   []                []; ...
                             'fileout'     'string'   []                EEG(1).filepath},...
                             'std_topo');
if ischar(g), error(g); end

% scan datasets
% -------------
all_topos = [];
for ind1 = 1:length(EEG)

    indMatch = [];
    for ind2 = 1:ind1-1
        if isequal(EEG(ind1).icawinv, EEG(ind2).icawinv)
            indMatch = ind2;
        end
    end
    
    if isempty(comps)
        compsToCompute = 1:size(EEG(ind1).icawinv,2);
    else
        compsToCompute = comps;
    end
    for k = compsToCompute
        if isempty(indMatch)
            % compute topo map grid (topoimage)
            % ---------------------------------
            chanlocs = EEG(ind1).chanlocs(EEG(ind1).icachansind);
            if isempty( [ chanlocs.theta ] )
                error('Channel locations are required for computing scalp topographies');
            end
            [~, grid, ~, Xi, Yi] = topoplot( EEG(ind1).icawinv(:,k), chanlocs, ...
                                                  'verbose', 'off',...
                                                   'electrodes', 'on' ,'style','both',...
                                                   'plotrad',0.55,'intrad',0.55,...
                                                   'noplot', 'on', 'chaninfo', EEG(ind1).chaninfo);
            all_topos = setfield(all_topos, [ 'comp' int2str(k) '_grid' ], { ':' ':' ind1 }, grid);
            all_topos = setfield(all_topos, [ 'comp' int2str(k) '_x' ]   , Xi(:,1));
            all_topos = setfield(all_topos, [ 'comp' int2str(k) '_y' ]   , Yi(:,1));
        else
            grid      = getfield(all_topos, [ 'comp' int2str(k) '_grid' ], { ':' ':' indMatch });
            all_topos = setfield(all_topos, [ 'comp' int2str(k) '_grid' ], { ':' ':' ind1 }, grid);
        end
    end
end

% Save topos in file
% ------------------
all_topos.datatype  = 'TOPO';
all_topos.trialinfo = g.trialinfo;
tmpfile = fullfile( g.fileout, [ EEG.filename(1:end-3) 'icatopo' ]); 
std_savedat(tmpfile, all_topos);
