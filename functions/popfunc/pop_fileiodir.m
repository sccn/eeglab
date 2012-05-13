% pop_fileiodir() - import directory into EEGLAB using FileIO 
%
% Usage:
%   >> OUTEEG = pop_fileiodir; % pop up window
%   >> OUTEEG = pop_fileiodir( folder );
%
% Inputs:
%   folder - [string] folder name
%
% Optional inputs:
%   'channels'   - [integer array] list of channel indices
%   'samples'    - [min max] sample point limits for importing data. 
%   'trials'     - [min max] trial's limit for importing data. 
%
% Outputs:
%   OUTEEG   - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2012-
%
% Note: FILEIO toolbox must be installed. 

% Copyright (C) 2012 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
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

function [EEG, command] = pop_fileiodir(folder, varargin); 
EEG = [];
command = '';

if nargin < 1
	% ask user
	folder = uigetdir('*.*', 'Choose a directory -- pop_fileiodir()'); 
	if folder == 0 return; end;
    drawnow;
    
    % open file to get infos
    % ----------------------
    disp('Reading data file header...');
    dat = ft_read_header(folder);
    uilist   = { { 'style' 'text' 'String' 'Channel list (defaut all):' } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'String' [ 'Data range (in sample points) (default all [1 ' int2str(dat.nSamples) '])' ] } ...
                 { 'style' 'edit' 'string' '' }  };
    geom = { [3 1] [3 1] };
    if dat.nTrials > 1
        uilist{end+1} = { 'style' 'text' 'String' [ 'Trial range (default all [1 ' int2str(dat.nTrials) '])' ] };
        uilist{end+1} = { 'style' 'edit' 'string' '' };
        geom = { geom{:} [3 1] };
    end;
    result = inputgui( geom, uilist, 'pophelp(''pop_fileiodir'')', 'Load data using FILE-IO -- pop_fileiodir()');
    if length(result) == 0 return; end;

    options = {};
    result = { result{:} '' };
    if ~isempty(result{1}), options = { options{:} 'channels' eval( [ '[' result{1} ']' ] ) }; end;
    if ~isempty(result{2}), options = { options{:} 'samples'  eval( [ '[' result{2} ']' ] ) }; end;
    if ~isempty(result{3}), options = { options{:} 'trials'   eval( [ '[' result{3} ']' ] ) }; end;
else
    dat = ft_read_header(folder);
    options = varargin;
end;

[EEG command] = pop_fileio(folder, options{:});

