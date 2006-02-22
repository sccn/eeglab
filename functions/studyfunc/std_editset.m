% editstudy() - modify study dataset.
%
% Usage: >> [STUDY, ALLEEG] = editstudy(STUDY, ALLEEG, key1, val1, ...);  
%
% Input:
%   STUDY      - STUDY set
%   ALLEEG     - EEGLAB vector of EEG sets included in the STUDY structure 
%
% Optional input:
%   'command'  - [cell array] change study (see command description and
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
%   'dipselect' - [float] select component with residual variance below
%                 the value given as input for all loaded datasets.
%
% Output:
%   STUDY      - a new STUDY set containing some or all of the datasets in ALLEEG, 
%                plus additional information from the optional inputs above. 
%   ALLEEG     - an EEGLAB vector of EEG sets included in the STUDY structure 
%
%  See also:  pop_createstudy(), load_alleeg(), pop_clust(), pop_preclust(), 
%             eeg_preclust(), eeg_createdata()
%
% Authors: Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, October , 2004-

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
% Revision 1.19  2006/02/09 21:48:13  arno
% implement rmclust
%
% Revision 1.18  2006/02/09 21:46:14  arno
% rmclust
%
% Revision 1.17  2006/02/09 19:44:26  arno
% same
%
% Revision 1.16  2006/02/09 19:43:45  arno
% scanning datasets
%
% Revision 1.15  2006/02/08 23:24:44  arno
% adding the comps entry
%
% Revision 1.14  2006/02/07 21:52:23  arno
% header
%
% Revision 1.13  2006/02/03 23:52:31  arno
% save field
%
% Revision 1.12  2006/02/03 23:16:08  arno
% adding .saved field
%
% Revision 1.11  2006/02/03 22:47:58  arno
% preserving save field when storing
%
% Revision 1.10  2006/02/03 22:33:41  arno
% fix updating the ALLEEG structure
%
% Revision 1.9  2006/02/03 22:12:50  arno
% fix typo
%
% Revision 1.8  2006/02/03 20:47:27  arno
% do not modify dat fields if they have not changed
%

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
                            'rmclust'   'string'  { 'on' 'off' }  'on';
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
            STUDY.datasetinfo(currentind).subject = g.commands{k+1};
        case 'comps'
            STUDY.datasetinfo(currentind).comps = g.commands{k+1};
        case 'condition'
            STUDY.datasetinfo(currentind).condition = g.commands{k+1};
        case 'group'
            STUDY.datasetinfo(currentind).group   = g.commands{k+1};
        case 'session' 
            STUDY.datasetinfo(currentind).session = g.commands{k+1};
        case 'remove'
            ALLEEG = eeg_store(ALLEEG, eeg_empty, g.commands{k+1});
        case 'dipselect'
            STUDY = checkstudy(STUDY); % update setind field
            
            for si = 1:size(STUDY.setind,2)% scan datasets that are part of STUDY
                
                % find a dataset with dipoles
                % ---------------------------
                idat = 0;
                for sc = 1:size(STUDY.setind,1)
                    if isfield(ALLEEG(STUDY.datasetinfo(STUDY.setind(sc,si)).index).dipfit, 'model')
                        idat = STUDY.datasetinfo(STUDY.setind(sc,si)).index;
                        break;
                    end;
                end;
                    
                if idat ~= 0
                    fprintf('Selecting dipole with less than %2.1f residual variance in dataset ''%s''\n', 100*rv, ALLEEG(idat).setname)
                    indleft = []; % components that are left in clustering
                    for icomp = succompind{si} % scan components
                        if (ALLEEG(idat).dipfit.model(icomp).rv < rv)
                             indleft = [indleft icomp];
                        end;
                    end;
                    STUDY.datasetinfo(STUDY.setind(sc,si)).comps = indleft;
                else
                    fprintf('No dipole information found in ''%s'' dataset, using all components\n', ALLEEG.setname)
                end
            end;
            
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
        if ALLEEG(currentind).session ~=          STUDY.datasetinfo(currentind).session
            ALLEEG(currentind).session          = STUDY.datasetinfo(currentind).session;
            ALLEEG(currentind).saved            = 'no';
        end;
        if ~strcmpi(char(ALLEEG(currentind).group), char(STUDY.datasetinfo(currentind).group))
            ALLEEG(currentind).group            = STUDY.datasetinfo(currentind).group;
            ALLEEG(currentind).saved            = 'no';
        end;
    end;
end;

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
[STUDY ALLEEG] = checkstudy(STUDY, ALLEEG);
if ~isempty(g.filename),
    [STUDY.filepath STUDY.filename ext] = fileparts(fullfile( g.filepath, g.filename ));
    STUDY.filename = [ STUDY.filename ext ];
    ver = version;
    STUDY.saved = 'yes';
    disp('Saving study...');
    if ver(1) > '6'
         save('-mat','-V6',fullfile( STUDY.filepath, STUDY.filename), 'STUDY');
    else save('-mat',      fullfile( STUDY.filepath, STUDY.filename), 'STUDY');
    end;
end
if strcmpi(g.resave, 'on')
    ver = version;
    disp('Saving study...');
    STUDY.saved = 'yes';
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
