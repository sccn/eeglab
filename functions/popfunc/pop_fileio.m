% pop_fileio() - import data files into EEGLAB using FileIO 
%
% Usage:
%   >> OUTEEG = pop_fileio; % pop up window
%   >> OUTEEG = pop_fileio( filename );
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'channels'   - [integer array] list of channel indices
%   'samples'    - [min max] sample point limits for importing data. 
%   'trials'     - [min max] trial's limit for importing data. 
%
% Outputs:
%   OUTEEG   - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2008-
%
% Note: FILEIO toolbox must be installed. 

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
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

function [EEG, command] = pop_fileio(filename, varargin); 
EEG = [];
command = '';

if nargin < 1
	% ask user
    ButtonName = questdlg2('Do you want to import a file or a folder?', ...
                           'FILE-IO import', ...
                           'Folder', 'File', 'File');
    if strcmpi(ButtonName, 'file')
        [filename, filepath] = uigetfile('*.*', 'Choose a file or header file -- pop_fileio()'); 
        drawnow;
        if filename(1) == 0 return; end;
        filename = fullfile(filepath, filename);
    else
        filename = uigetfile('*.*', 'Choose a folder -- pop_fileio()'); 
        drawnow;
        if filename(1) == 0 return; end;
    end;
    
    % open file to get infos
    % ----------------------
    disp('Reading data file header...');
    dat = ft_read_header(filename);
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
    result = inputgui( geom, uilist, 'pophelp(''pop_fileio'')', 'Load data using FILE-IO -- pop_fileio()');
    if length(result) == 0 return; end;

    options = {};
    result = { result{:} '' };
    if ~isempty(result{1}), options = { options{:} 'channels' eval( [ '[' result{1} ']' ] ) }; end;
    if ~isempty(result{2}), options = { options{:} 'samples'  eval( [ '[' result{2} ']' ] ) }; end;
    if ~isempty(result{3}), options = { options{:} 'trials'   eval( [ '[' result{3} ']' ] ) }; end;
else
    dat = ft_read_header(filename);
    options = varargin;
end;

% decode imput parameters
% -----------------------
g = finputcheck( options, { 'samples'     'integer' [1 Inf]    [];
                            'trials'      'integer' [1 Inf]    [];
                            'channels'    'integer' [1 Inf]    [] }, 'pop_fileio');
if isstr(g), error(g); end;

% import data
% -----------
EEG = eeg_emptyset;
fprintf('Reading data ...\n');
dataopts = {};
if ~isempty(g.channels), dataopts = { dataopts{:} 'chanindx', g.channels }; end;
if ~isempty(g.samples ), dataopts = { dataopts{:} 'begsample', g.samples(1), 'endsample', g.samples(2)}; end;
if ~isempty(g.trials  ), dataopts = { dataopts{:} 'begtrial', g.trials(1), 'endtrial', g.trials(2)}; end;
alldata = ft_read_data(filename, 'header', dat, dataopts{:});

% convert to seconds for sread
% ----------------------------
EEG.srate           = dat.Fs;
EEG.nbchan          = dat.nChans;
EEG.data            = alldata;
EEG.setname 		= '';
EEG.comments        = [ 'Original file: ' filename ];
EEG.xmin = -dat.nSamplesPre/EEG.srate; 
EEG.trials          = dat.nTrials;
EEG.pnts            = dat.nSamples;
if isfield(dat, 'label') && ~isempty(dat.label)
    EEG.chanlocs = struct('labels', dat.label);
end

% extract events
% --------------
disp('Reading events...');
try
    event = ft_read_event(filename);
catch, disp(lasterr); event = []; end;
if ~isempty(event)
    subsample = 0;
    
    if ~isempty(g.samples), subsample = g.samples(1); end;
    
    for index = 1:length(event)
        offset = fastif(isempty(event(index).offset), 0, event(index).offset);
        EEG.event(index).type     = event(index).value;
        EEG.event(index).value    = event(index).type;
        EEG.event(index).latency  = event(index).sample+offset+subsample;
        EEG.event(index).duration = event(index).duration;
        if EEG.trials > 1
            EEG.event(index).epoch = ceil(EEG.event(index).latency/EEG.pnts);        
        end;
    end;
    
    EEG = eeg_checkset(EEG, 'eventconsistency');
else 
    disp('Warning: no event found. Events might be embeded in a data channel.');
    disp('         To extract events, use menu File > Import Event Info > From data channel');
end;

% convert data to single if necessary
% -----------------------------------
EEG = eeg_checkset(EEG,'makeur');   % Make EEG.urevent field

% history
% -------
if isempty(options)
    command = sprintf('EEG = pop_fileio(''%s'');', filename); 
else
    command = sprintf('EEG = pop_fileio(''%s'', %s);', filename, vararg2str(options)); 
end;    
