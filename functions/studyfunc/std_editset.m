% editstudy() - modify study dataset.
%
% Usage: >> [STUDY, ALLEEG] = editstudy(STUDY, ALLEEG, key1, val1, ...);  
%
% Input:
%   STUDY      - STUDY set
%   ALLEEG     - EEGLAB vector of EEG sets included in the STUDY structure 
%
% Optional input:
%   'command'   - [cell array] change study (see command description and
%                 example below.
%   'name'     - [string] a specified (mnemonic) name for the STUDY structure. 
%                {default: ''}
%   'task'     - [string] a description of the experimental task(s) performed 
%                by the STUDY subjects {default: ''}.  
%   'filename' - [string] filename for the STUDY set.
%   'filepath' - [string] file path (directory/folder) in which the STUDY file
%                will be saved.  
%   'notes'    - [string] notes about the experiment, the datasets, the STUDY, 
%                or anything to keep in the record {default: ''}. 
%   'updatedat' - ['on'|'off'] update 'subject' 'session' 'condition' and 
%                group field of datasets.
%   'savedat'   - ['on'|'off'] resave datasets
%
% Each command is a cell array composed of the following: 
%   'index'     - [integer] mmodify dataset index.
%   'remove'    - [integer] remove dataset index.
%   'subject'   - [string] subject name or code.
%   'condition' - [string] condition corresponding to dataset.
%   'session '  - [integer] dataset session index.
%   'group'     - [string] dataset group.
%   'subject'   - [string] subject name or code.
%   'load'      - [string] load dataset having filename contained in entry.
%
% Output:
%   STUDY      - a new STUDY set containing some or all of the datasets in ALLEEG, 
%                plus additional information from the optional inputs above. 
%   ALLEEG     - an EEGLAB vector of EEG sets included in the STUDY structure 
%
%  See also:  pop_createstudy(), load_ALLEEG(), pop_clust(), pop_preclust(), 
%             eeg_preclust(), eeg_createdata()
%
% Authors:  Hilit Serby, Arnaud Delorme, SCCN, INC, UCSD, October , 2004-

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

function [STUDY, ALLEEG] = editstudy(STUDY, ALLEEG, varargin) 

if (nargin < 3)
    help editstudy;
    return;
end;

% decode input parameters
% -----------------------
g = finputcheck(varargin, { 'updatedat' 'string'  { 'on' 'off' }  'on';
                            'name'      'string'  { }             '';
                            'task'      'string'  { }             '';
                            'notes'     'string'  { }             '';
                            'filename'  'string'  { }             '';
                            'filepath'  'string'  { }             '';
                            'resave'    'string'  { 'on' 'off' 'info' }  'off';
                            'savedat'   'string'  { 'on' 'off' }  'on';
                            'commands'  'cell'    {}              {} }, 'editstudy');
if isstr(g), error(g); end;

if ~isempty(g.name),  STUDY.name  = g.name; end
if ~isempty(g.task),  STUDY.task  = g.task; end
if ~isempty(g.notes), STUDY.notes = g.notes; end

% make one cell array with commands
% ---------------------------------
allcoms = {};
if ~isempty(g.commands)
    if iscell(g.commands{1})
        for k = 1:length(g.commands)
            allcoms = { allcoms{:} g.commands{k}{:} };
        end;
    else 
        allcoms = g.commands;
    end;
end;
g.commands = allcoms;

% copy values
% -----------
if ~isfield(STUDY, 'datasetinfo')
    for realindex = 1:length(ALLEEG)
        if ~isempty(ALLEEG(realindex).data)
            [tmppath tmpfile tmpext] = fileparts(  fullfile(ALLEEG(realindex).filepath, ALLEEG(realindex).filename) );
            STUDY.datasetinfo(realindex).filepath  = tmppath;   
            STUDY.datasetinfo(realindex).filename  = [ tmpfile tmpext ];   
            STUDY.datasetinfo(realindex).subject   = ALLEEG(realindex).subject;
            STUDY.datasetinfo(realindex).session   = ALLEEG(realindex).session;
            STUDY.datasetinfo(realindex).condition = ALLEEG(realindex).condition;
            STUDY.datasetinfo(realindex).group     = ALLEEG(realindex).group;                    
        end;
    end;
end;
        
% execute commands
% ----------------
currentind = 1;
rmlist = [];
for k = 1:2:length(g.commands)
    switch g.commands{k}
        case 'index'
            currentind = g.commands{k+1};
        case 'subject'
            if strcmpi(g.updatedat, 'on')
                if ~strcmpi(ALLEEG(currentind).subject, g.commands{k+1})
                    ALLEEG(currentind).subject        = g.commands{k+1};
                    ALLEEG(currentind).saved          = 'no';
                end;
            end; 
            STUDY.datasetinfo(currentind).subject = g.commands{k+1};
        case 'condition'
            if strcmpi(g.updatedat, 'on')
                if ~strcmpi(ALLEEG(currentind).condition, g.commands{k+1})
                    ALLEEG(currentind).condition      = g.commands{k+1};
                    ALLEEG(currentind).saved          = 'no';
                end;
            end; 
            STUDY.datasetinfo(currentind).condition = g.commands{k+1};
        case 'group'
            if strcmpi(g.updatedat, 'on')
                if ~strcmpi(ALLEEG(currentind).group, g.commands{k+1})
                    ALLEEG(currentind).group          = g.commands{k+1};
                    ALLEEG(currentind).saved          = 'no';
                end;
            end; 
            STUDY.datasetinfo(currentind).group   = g.commands{k+1};
        case 'session' 
            if strcmpi(g.updatedat, 'on')
                if session(ALLEEG(currentind).session ~= g.commands{k+1}
                    ALLEEG(currentind).session        = g.commands{k+1};
                    ALLEEG(currentind).saved          = 'no';
                end;
            end; 
            STUDY.datasetinfo(currentind).session = g.commands{k+1};
        case 'remove'
            ALLEEG = eeg_store(ALLEEG, eeg_empty, g.commands{k+1});
        case 'load'
            TMPEEG = pop_loadset('filename', g.commands{k+1}, 'loadmode', 'info');
            ALLEEG = eeg_store(ALLEEG, eeg_checkset(TMPEEG), currentind);
            
            % update datasetinfo structure
            % ----------------------------
            [tmppath tmpfile tmpext] = fileparts( fullfile(ALLEEG(currentind).filepath, ...
                                                           ALLEEG(currentind).filename) );
            STUDY.datasetinfo(currentind).filepath  = tmppath;   
            STUDY.datasetinfo(currentind).filename  = [ tmpfile tmpext ];   
            STUDY.datasetinfo(currentind).subject   = ALLEEG(currentind).subject;
            STUDY.datasetinfo(currentind).session   = ALLEEG(currentind).session;
            STUDY.datasetinfo(currentind).condition = ALLEEG(currentind).condition;
            STUDY.datasetinfo(currentind).group     = ALLEEG(currentind).group;                    
    end
end

% remove empty datasets (cnnot be done above because some empty datasets might not have been removed
% ---------------------
[ ALLEEG STUDY.datasetinfo ] = removeempty(ALLEEG, STUDY.datasetinfo);

% save datasets if necessary
% --------------------------
if strcmpi(g.savedat, 'on')
    for index = 1:length(ALLEEG)
        if isempty(ALLEEG(index).filename)
            fprintf('Cannot resave ALLEEG(%d) because the dataset has no filename\n', index);
        else
            TMP = pop_saveset(ALLEEG(index), 'savemode', 'resave');
            ALLEEG = eeg_store(ALLEEG, TMP, index);
        end;
    end;
end;

% save study if necessary
% -----------------------
[STUDY ALLEEG] = checkstudy(STUDY, ALLEEG);
if ~isempty(g.filename),
    [STUDY.filepath STUDY.filename ext] = fileparts(fullfile( g.filepath, g.filename ));
    STUDY.filename = [ STUDY.filename ext ];
    ver = version;
    disp('Saving study...');
    if ver(1) > '6'
         save('-mat','-V6',fullfile( STUDY.filepath, STUDY.filename), 'STUDY');
    else save('-mat',      fullfile( STUDY.filepath, STUDY.filename), 'STUDY');
    end;
end
if strcmpi(g.resave, 'on')
    ver = version;
    disp('Saving study...');
    if ver(1) > '6'
         save('-mat','-V6',fullfile( STUDY.filepath, STUDY.filename), 'STUDY');
    else save('-mat',      fullfile( STUDY.filepath, STUDY.filename), 'STUDY');
    end;
end;    

% ---------------------
% remove empty elements
% ---------------------
function [ALLEEG, datasetinfo] = removeempty(ALLEEG, datasetinfo);

    rmindex = [];
    for index = 1:length(datasetinfo)
        if isempty(datasetinfo(index).subject) & isempty(ALLEEG(index).nbchan)
            rmindex = [ rmindex index ];
        end;
    end;
    datasetinfo(rmindex) = [];
    ALLEEG(rmindex)      = [];
