% pop_loadstudy() - load an existing EEGLAB STUDY set of EEG datasets plus 
%                   its corresponding ALLEEG structure. Calls std_loadalleeg().
% Usage:
%   >> [STUDY ALLEEG] = pop_loadstudy; % pop up a window to collect filename
%   >> [STUDY ALLEEG] = pop_loadstudy( 'key', 'val', ...); % no pop-up
%
% Optional inputs:
%   'filename' - [string] filename of the STUDY set file to load.
%   'filepath' - [string] filepath of the STUDY set file to load.
%
% Outputs:
%   STUDY      - the requested STUDY set structure.
%   ALLEEG     - the corresponding ALLEEG structure containing 
%                the (loaded) STUDY EEG datasets.    
%
% See also: std_loadalleeg(), pop_savestudy()
%
% Authors: Hilit Serby & Arnaud Delorme, SCCN, INC, UCSD, September 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, Spetember 2005, hilit@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.


function [STUDY, ALLEEG, com] = pop_loadstudy(varargin)

STUDY  = [];
ALLEEG = [];
com = '';
if isempty(varargin)
    [filename, filepath] = uigetfile2('*.study', 'Load a STUDY -- pop_loadstudy()'); 
    if filename(1) == 0, return; end;
    if ~strncmp(filename(end-5:end), '.study',6)
        if isempty(strfind(filename,'.'))
            filename = [filename '.study'];
        else
            filename = [filename(1:strfind(filename,'.')-1) '.study'];
        end
    end
else
    filepath = '';
    if nargin == 1
        varargin = { 'filename' varargin{:} };
    end;
    for k = 1:2:length(varargin)
        switch varargin{k}
            case 'filename'
                filename = varargin{k+1};
            case 'filepath'
                filepath = varargin{k+1};
        end
    end
end

if ~isempty(filename)
    STUDYfile = fullfile(filepath,filename);
    try 
        load('-mat', STUDYfile);
    catch
        error(['pop_loadstudy(): STUDY set file ''STUDYfile'' not loaded -- check filename and path']);
    end
    [filepath filename ext] = fileparts(STUDYfile);
    STUDY.filename = [filename ext];
    STUDY.etc.oldfilepath = STUDY.filepath;
    STUDY.filepath        = filepath;
else
    error(['pop_loadstudy(): No STUDY set file provided.']);
end
  
ALLEEG = std_loadalleeg(STUDY);

% Update the pointers from STUDY to the ALLEEG datasets
for k = 1:length(STUDY.datasetinfo)
    STUDY.datasetinfo(k).index = k;
    STUDY.datasetinfo(k).filename = ALLEEG(k).filename;
    STUDY.datasetinfo(k).filepath = ALLEEG(k).filepath;
end

if ~isfield(STUDY, 'changrp'), STUDY.changrp = []; end;
[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);

if ~isfield(STUDY, 'changrp') || isempty(STUDY.changrp)
    if std_uniformfiles(STUDY, ALLEEG) == 0
         STUDY = std_changroup(STUDY, ALLEEG);
    else STUDY = std_changroup(STUDY, ALLEEG, [], 'interp');
    end;
end;

% Update the design path
for inddes = 1:length(STUDY.design)
    for indcell = 1:length(STUDY.design(inddes).cell)
        if isempty(STUDY.design(inddes).filepath)
            pathname = STUDY.datasetinfo(STUDY.design(inddes).cell(indcell).dataset(1)).filepath;
        else
            pathname = STUDY.design(inddes).filepath;
        end
        filebase = STUDY.design(inddes).cell(indcell).filebase;
        tmpinds1 = find(filebase == '/');
        tmpinds2 = find(filebase == '\');
        if ~isempty(tmpinds1)
            STUDY.design(inddes).cell(indcell).filebase = fullfile(pathname, filebase(tmpinds1(end)+1:end));
        elseif ~isempty(tmpinds2)
            STUDY.design(inddes).cell(indcell).filebase = fullfile(pathname, filebase(tmpinds2(end)+1:end));
        else STUDY.design(inddes).cell(indcell).filebase = fullfile(pathname, filebase );
        end;
    end;
end;

% check for corrupted ERSP ICA data files
% A corrupted file is present if
% - components have been selected
% - .icaersp or .icaitc files are present
% - the .trialindices field is missing from these files
try
    %% check for corrupted ERSP ICA data files
    ncomps1 = cellfun(@length, { STUDY.datasetinfo.comps });
    ncomps2 = cellfun(@(x)(size(x,1)), { ALLEEG.icaweights });
    if any(~isempty(ncomps1))
        if any(ncomps1 ~= ncomps2)
            warningshown = 0;

            for des = 1:length(STUDY.design)
                for iCell = 1:length(STUDY.design(des).cell)
                    if ~warningshown
                        if exist( [ STUDY.design(des).cell(iCell).filebase '.icaersp' ] )
                            warning('off', 'MATLAB:load:variableNotFound');
                            tmp = load('-mat', [ STUDY.design(des).cell(iCell).filebase '.icaersp' ], 'trialindices');
                            warning('on', 'MATLAB:load:variableNotFound');
                            if ~isfield(tmp, 'trialindices')
                                warningshown = 1;
                                warndlg( [ 'Warning: ICA ERSP or ITC data files computed with old version of EEGLAB for design ' int2str(des) 10 ...
                                             '(and maybe other designs). These files may be corrupted and must be recomputed.' ], 'Important EEGLAB warning', 'nonmodal');
                            end;
                        end;
                        if warningshown == 0 && exist( [ STUDY.design(des).cell(iCell).filebase '.icaitc' ] )
                            tmp = load('-mat', [ STUDY.design(des).cell(iCell).filebase '.icaersp' ], 'trialindices');
                            if ~isfield(tmp, 'trialindices')
                                warningshown = 1;
                                warndlg( [ 'Warning: ICA ERSP or ITC data files computed with old version of EEGLAB for design ' int2str(des) 10 ...
                                             '(and maybe other designs). These files may be corrupted and must be recomputed.' ], 'Important EEGLAB warning', 'modal');
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
catch, 
    disp('Warning: failed to test STUDY file version');
end;

TMP = STUDY.datasetinfo;
STUDY = std_maketrialinfo(STUDY, ALLEEG);
if ~isequal(STUDY.datasetinfo, TMP)
    disp('STUDY Warning: the trial information collected from datasets has changed');
end;
std_checkfiles(STUDY, ALLEEG);
STUDY.saved = 'yes';
STUDY = std_selectdesign(STUDY, ALLEEG, STUDY.currentdesign);

com = sprintf('[STUDY ALLEEG] = pop_loadstudy(''filename'', ''%s'', ''filepath'', ''%s'');', STUDY.filename, STUDY.filepath);
