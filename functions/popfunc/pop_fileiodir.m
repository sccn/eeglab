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

function [EEG, command] = pop_fileiodir(folder, varargin); 
EEG = [];
command = '';

if nargin < 1
	% ask user
	folder = uigetdir('*.*', 'Choose a directory -- pop_fileiodir()'); 
	if folder == 0 return; end
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
    end
    result = inputgui( geom, uilist, 'pophelp(''pop_fileiodir'')', 'Load data using FILE-IO -- pop_fileiodir()');
    if length(result) == 0 return; end

    options = {};
    result = { result{:} '' };
    if ~isempty(result{1}), options = { options{:} 'channels' eval( [ '[' result{1} ']' ] ) }; end
    if ~isempty(result{2}), options = { options{:} 'samples'  eval( [ '[' result{2} ']' ] ) }; end
    if ~isempty(result{3}), options = { options{:} 'trials'   eval( [ '[' result{3} ']' ] ) }; end
else
    dat = ft_read_header(folder);
    options = varargin;
end

[EEG command] = pop_fileio(folder, options{:});

