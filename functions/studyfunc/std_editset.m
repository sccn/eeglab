% std_editset() - modify a STUDY set structure.
%
% Usage: 
%             >> [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, key1, val1, ...);  
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   ALLEEG     - vector of the EEG datasets included in the STUDY structure 
%
% Optional inputs:
%   'commands' - {cell_array} change STUDY (see command description and
%                 example below.
%   'name'     - [string] specify a (mnemonic) name for the STUDY structure. 
%                {default: ''}
%   'task'     - [string] attach a description of the experimental task(s) 
%                performed by the STUDY subjects {default: ''}.  
%   'filename' - [string] filename for the STUDY set.
%   'filepath' - [string] file path (directory/folder) in which the STUDY file
%                will be saved.  
%   'addchannellabels' - ['on'|'off'] add channel labels ('1', '2', '3', ...)
%                to all datasets of a STUDY to ensure that all STUDY functions
%                will work {default: 'off' unless no dataset has channel
%                locations and then it is automatically set to on}
%   'notes'    - [string] notes about the experiment, the datasets, the STUDY, 
%                or anything else to store with the STUDY itself {default: ''}. 
%   'updatedat' - ['on'|'off'] update 'subject' 'session' 'condition' and/or
%                'group' fields of STUDY dataset(s).
%   'savedat'   - ['on'|'off'] re-save datasets
%   'inbrain'   - ['on'|'off'] select components for clustering from all STUDY 
%                 datasets with equivalent dipoles located inside the brain volume. 
%                 Dipoles are selected based on their residual variance and their 
%                 location {default: 'off'}
%   'resave'    - ['on'|'off'] save or resave STUDY {default: 'off'}
%
% Each of the 'commands' (above) is a cell array composed of any of the following: 
%   'index'     - [integer] modify/add dataset index. Note that if a
%                 dataset is added and that this leaves some indices not 
%                 populated, the dataset is automatically set to the last
%                 empty index. For instance creating a STUDY with a single
%                 dataset at index 10 will result with a STUDY with a
%                 single dataset at index 1.
%   'remove'    - [integer] remove dataset index.
%   'subject'   - [string] subject code.
%   'condition' - [string] dataset condition. 
%   'session '  - [integer] dataset session number.
%   'group'     - [string] dataset group.
%   'load'      - [filename] load dataset from specified filename 
%   'dipselect' - [float<1] select components for clustering from all STUDY 
%                 datasets with dipole model residual var. below this value. 
%   'inbrain'   - ['on'|'off'] same as above. This option may also be
%                 placed in the command list (preceeding the 'dipselect'
%                 option).
%
% Outputs:
%   STUDY      - a new STUDY set containing some or all of the datasets in ALLEEG, 
%                plus additional information from the optional inputs above. 
%   ALLEEG     - a vector of EEG datasets included in the STUDY structure 
%
%  See also:  pop_createstudy(), std_loadalleeg(), pop_clust(), pop_preclust(), 
%             eeg_preclust(), eeg_createdata()
%
% Authors: Arnaud Delorme, Hilit Serby, SCCN/INC/UCSD, October , 2004-

% Copyright (C) Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD, October 11, 2004, smakeig@ucsd.edu
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

function [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, varargin) 

if (nargin < 3)
    help std_editset;
    return;
end;

% decode input parameters
% -----------------------
g = finputcheck(varargin, { 'updatedat' 'string'  { 'on','off' }  'off';
                            'name'      'string'  { }             '';
                            'task'      'string'  { }             '';
                            'notes'     'string'  { }             '';
                            'filename'  'string'  { }             '';
                            'filepath'  'string'  { }             '';
                            'resave'    'string'  { 'on','off','info' }  'off';
                            'savedat'   'string'  { 'on','off' }  'off';
                            'addchannellabels' 'string'  { 'on','off' }  'off';
                            'rmclust'   'string'  { 'on','off' }  'on';
                            'inbrain'   'string'  { 'on','off' }  'off';
                            'commands'  'cell'    {}              {} }, 'std_editset');
if isstr(g), error(g); end;

if isempty(STUDY), STUDY.history = 'STUDY = [];'; end;
if ~isempty(g.name),  STUDY.name  = g.name; end
if ~isempty(g.task),  STUDY.task  = g.task; end
if ~isempty(g.notes), STUDY.notes = g.notes; end

% default addchannellabels
% ------------------------
if ~isempty(ALLEEG)
    allchanlocs = { ALLEEG.chanlocs };
    if all(cellfun( @isempty, allchanlocs))
        g.addchannellabels = 'on';
    else
        if any(cellfun( @isempty, allchanlocs))
            error( [ 'Some datasets have channel locations and some other don''t' 10 ...
                     'the STUDY is not homogenous and cannot be created.' ]);
        end;
    end;
end;

% make one cell array with commands
% ---------------------------------
allcoms = {};
if ~isempty(g.commands)
    if iscell(g.commands{1})
        for k = 1:length(g.commands)
            % put index field first
            indindex   = strmatch('index', lower(g.commands{k}(1:2:end)));
            if ~isempty(indindex)
                 tmpcom = { 'index' g.commands{k}{2*(indindex-1)+1+1} g.commands{k}{:} };
            else tmpcom = g.commands{k}; 
            end;
            allcoms = { allcoms{:} tmpcom{:} };
        end;
    else 
        allcoms = g.commands;
    end;
end;
g.commands = allcoms;

% add 'dipselect' command if 'inbrain' option is selected
% ---------------------------------
dipselectExists = false;
for k = 1:2:length(g.commands)
    if strcmp(g.commands{k},'dipselect')
       dipselectExists = true;
    end;
end;
if strcmp(g.inbrain,'on') && ~dipselectExists
    g.commands{length(g.commands)+1} = 'dipselect';
    g.commands{length(g.commands)+1} = 0.15;
end;

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
            STUDY.datasetinfo(currentind).subject = g.commands{k+1};
        case 'comps'
            STUDY.datasetinfo(currentind).comps = g.commands{k+1};
        case 'condition'
            STUDY.datasetinfo(currentind).condition = g.commands{k+1};
        case 'group'
            STUDY.datasetinfo(currentind).group   = g.commands{k+1};
        case 'session' 
            STUDY.datasetinfo(currentind).session = g.commands{k+1};
        case 'session' 
            STUDY.datasetinfo(currentind).session = g.commands{k+1};
        case 'remove'
            % create empty structure
            allfields = fieldnames(ALLEEG);
            tmpfields = allfields;
            tmpfields(:,2) = cell(size(tmpfields));
            tmpfields = tmpfields';
            ALLEEG(g.commands{k+1}) = struct(tmpfields{:});

            % create empty structure
            allfields = fieldnames(STUDY.datasetinfo);
            tmpfields = allfields;
            tmpfields(:,2) = cell(size(tmpfields));
            tmpfields = tmpfields';
            STUDY.datasetinfo(g.commands{k+1}) = struct(tmpfields{:});
            
            if isfield(STUDY.datasetinfo, 'index')
                STUDY.datasetinfo = rmfield(STUDY.datasetinfo, 'index');
            end;
            STUDY.datasetinfo(1).index = [];
            STUDY.changrp = [];
        case 'return', return;
        case 'inbrain' 
            g.inbrain = g.commands{k+1};
        case 'dipselect'
            STUDY = std_checkset(STUDY, ALLEEG);
            rv = g.commands{k+1};
            clusters = std_findsameica(ALLEEG);
            
            for cc = 1:length(clusters)
                
                idat = 0;
                for tmpi = 1:length(clusters{cc})
                    if isfield(ALLEEG(clusters{cc}(tmpi)).dipfit, 'model')
                        idat = clusters{cc}(tmpi);
                    end;
                end;
                
                indleft = [];
                if rv ~= 1
                    if idat ~= 0
                        if strcmp(g.inbrain,'on')
                            fprintf('Selecting dipoles with less than %%%2.1f residual variance and removing dipoles outside brain volume in dataset ''%s''\n', ...
                                100*rv, ALLEEG(idat).setname);
                            indleft = eeg_dipselect(ALLEEG(idat), rv*100,'inbrain'); 
                        else
                           fprintf('Selecting dipoles with less than %%%2.1f residual variance in dataset ''%s''\n', ...
                                100*rv, ALLEEG(idat).setname);
                            indleft = eeg_dipselect(ALLEEG(idat), rv*100,'rv'); 
                        end;
                    else
                        fprintf('No dipole information found in ''%s'' dataset, using all components\n', ALLEEG.setname)
                    end
                end;
                for tmpi = 1:length(clusters{cc})                
                    STUDY.datasetinfo(clusters{cc}(tmpi)).comps = indleft;
                end;
            end;
            STUDY.cluster = [];
            STUDY = std_checkset(STUDY, ALLEEG); % recreate parent dataset           
           
        case 'load'
            TMPEEG = std_loadalleeg( { g.commands{k+1} } );
            ALLEEG = eeg_store(ALLEEG, eeg_checkset(TMPEEG), currentind);
            ALLEEG(currentind).saved = 'yes';
            
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
            STUDY.datasetinfo(currentind).index     = currentind;    
        otherwise, error(sprintf('Unknown command %s', g.commands{k}));
    end
end

% add channel labels automatically
% -------------------------------
if strcmpi(g.addchannellabels, 'on')
    disp('Generating channel labels for all datasets...');
    for currentind = 1:length(ALLEEG)
        for ind = 1:ALLEEG(currentind).nbchan
            ALLEEG(currentind).chanlocs(ind).labels = int2str(ind);
        end;
    end;
    ALLEEG(currentind).saved = 'no';
    g.savedat = 'on';
end;
    
% update ALLEEG structure?
% ------------------------
if strcmpi(g.updatedat, 'on')
    for currentind = 1:length(ALLEEG)
        if ~strcmpi(ALLEEG(currentind).subject,   STUDY.datasetinfo(currentind).subject)
            ALLEEG(currentind).subject          = STUDY.datasetinfo(currentind).subject;
            ALLEEG(currentind).saved            = 'no';
        end;
        if ~strcmpi(ALLEEG(currentind).condition, STUDY.datasetinfo(currentind).condition)
            ALLEEG(currentind).condition        = STUDY.datasetinfo(currentind).condition;
            ALLEEG(currentind).saved            = 'no';
        end;
        if ~isequal(ALLEEG(currentind).session, STUDY.datasetinfo(currentind).session)
            ALLEEG(currentind).session          = STUDY.datasetinfo(currentind).session;
            ALLEEG(currentind).saved            = 'no';
        end;
        if ~strcmpi(char(ALLEEG(currentind).group), char(STUDY.datasetinfo(currentind).group))
            ALLEEG(currentind).group            = STUDY.datasetinfo(currentind).group;
            ALLEEG(currentind).saved            = 'no';
        end;
    end;
end;

% remove empty datasets (cannot be done above because some empty datasets
% might not have been removed)
% ---------------------
rmindex = [];
for index = 1:length(STUDY.datasetinfo)
    if isempty(STUDY.datasetinfo(index).subject) && isempty(ALLEEG(index).nbchan)
        rmindex = [ rmindex index ];
    end;
end;
STUDY.datasetinfo(rmindex) = [];
ALLEEG(rmindex)            = [];
for index = 1:length(STUDY.datasetinfo)
    STUDY.datasetinfo(index).index = index;
end;

% remove empty ALLEEG structures
% ------------------------------
while length(ALLEEG) > length(STUDY.datasetinfo)
   ALLEEG(end) = [];
end;
%[ ALLEEG STUDY.datasetinfo ] = remove_empty(ALLEEG, STUDY.datasetinfo);

% save datasets if necessary
% --------------------------
if strcmpi(g.savedat, 'on')
    for index = 1:length(ALLEEG)
        if isempty(ALLEEG(index).filename)
            fprintf('Cannot resave ALLEEG(%d) because the dataset has no filename\n', index);
        else
            TMP = pop_saveset(ALLEEG(index), 'savemode', 'resave');
            ALLEEG = eeg_store(ALLEEG, TMP, index);
            ALLEEG(index).saved = 'yes';
        end;
    end;
end;

% remove cluster information if necessary
% ---------------------------------------
if strcmpi(g.rmclust, 'on')
    STUDY.cluster = [];
end;

% save study if necessary
% -----------------------
if ~isempty(g.commands)
    STUDY.changrp = [];
    STUDY.cluster = [];
%    if ~isempty(STUDY.design)
%        [STUDY] = std_createclust(STUDY, ALLEEG, 'parentcluster', 'on');
%    end;
end;
[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);
if ~isempty(g.filename),
    [STUDY.filepath STUDY.filename ext] = fileparts(fullfile( g.filepath, g.filename ));
    STUDY.filename = [ STUDY.filename ext ];
    g.resave = 'on';
end
if strcmpi(g.resave, 'on')
    STUDY = pop_savestudy(STUDY, ALLEEG, 'savemode', 'resave');
end;    

% ---------------------
% remove empty elements
% ---------------------
function [ALLEEG, datasetinfo] = remove_empty(ALLEEG, datasetinfo);

    rmindex = [];
    for index = 1:length(datasetinfo)
        if isempty(datasetinfo(index).subject) && isempty(ALLEEG(index).nbchan)
            rmindex = [ rmindex index ];
        end;
    end;
    datasetinfo(rmindex) = [];
    ALLEEG(rmindex)      = [];
    for index = 1:length(datasetinfo)
        datasetinfo(index).index = index;
    end;
    
    % remove empty ALLEEG structures
    % ------------------------------
    while length(ALLEEG) > length(datasetinfo)
       ALLEEG(end) = [];
    end;
        
    
    
    
