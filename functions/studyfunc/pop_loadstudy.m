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

% Coding notes: Useful information on functions and global variables used.


function [STUDY, ALLEEG, com] = pop_loadstudy(varargin)

STUDY  = [];
ALLEEG = [];
com = '';
if isempty(varargin)
    [filename, filepath] = uigetfile2('*.study', 'Load a STUDY -- pop_loadstudy()'); 
    if filename(1) == 0, return; end
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
    end
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

% check for old study format
if ~isfield(STUDY, 'design')
    STUDY.design = [];
end
if ~isempty(STUDY.design)
    if isfield(STUDY.design, 'cell')
        txt = [ 'You are loading a STUDY from a previous version of EEGLAB.' 10 ...
                'Most study settings are backward compatible but require that' 10 ... 
                'you recompute the measure data files for all measures (ERP, ERSP,' 10 ...
                'ITC, spectrum).' ];
        dbs = dbstack;
        if isempty(varargin) % means that it was called from a call back
            warndlg2(txt);
        else
            fprintf(2,[txt 10]);
        end
    end
end

% check if old study format has different subjects
if ~isempty(STUDY.design)
    cases = {STUDY.design.cases};
    if ~all(cellfun(@(x)isequal(x, cases{end}), cases))
        dbs = dbstack;
        txt = [ 'You are loading a STUDY from a previous version of EEGLAB' 10 ...
                'that has designs with different subjects included. It is' 10 ....
                'no longer possible to have different subjects included in' 10 ...
                'different designs (to do so create a separate study). All' 10 ....
                'designs have been changed to include the same subjects' 10 ...
                'as design number 1.' ];
        if isempty(varargin) == 1 % means that it was called from a call back
            warndlg2(txt);
        else
            fprintf(2,[txt 10]);
        end
        for iDes = 2:length(STUDY.design)
            STUDY.design(iDes).cases = STUDY.design(1).cases;
        end
    end
end                      

if ~isfield(STUDY, 'changrp'), STUDY.changrp = []; end
if isempty(varargin)
     [STUDY ALLEEG] = std_checkset(STUDY, ALLEEG, 'popup');
else [STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);
end

if ~isfield(STUDY, 'changrp') || isempty(STUDY.changrp)
    if std_uniformfiles(STUDY, ALLEEG) == 0
         STUDY = std_changroup(STUDY, ALLEEG);
    else STUDY = std_changroup(STUDY, ALLEEG, [], 'interp');
    end
end

TMP = STUDY.datasetinfo;
STUDYTMP = std_maketrialinfo(STUDY, ALLEEG); % some dataset do not have trialinfo and
sameTrialInfo = true;
if isfield(STUDY.datasetinfo, 'trialinfo') ~= isfield(STUDYTMP.datasetinfo, 'trialinfo')
    sameTrialInfo = false;
elseif isfield(STUDY.datasetinfo, 'trialinfo') && isfield(STUDYTMP.datasetinfo, 'trialinfo')
    if exist('isequaln', 'builtin')
        sameTrialInfo = isequaln( { STUDY.datasetinfo.trialinfo }, { STUDYTMP.datasetinfo.trialinfo });
    else
        sameTrialInfo = true;
    end
end
if ~sameTrialInfo
    if isempty(varargin)
        res = questdlg2([ 'Warning: trial information contained in datasets not consistent' 10 'with STUDY information. Should EEGLAB overwrite STUDY trial information?' ], 'STUDY Warning', 'No', 'Yes', 'Yes');
        if isequal(res, 'Yes')
            STUDY = STUDYTMP;
        end
    else
        disp('Warning: trial information contained in datasets not consistent with STUDY information.');
        disp('         In rare cases, this could be intentional. Use "STUDY = std_maketrialinfo(STUDY, ALLEEG);" to fix.');
    end
end
clear STUDYTMP;
std_checkfiles(STUDY, ALLEEG);
STUDY.saved = 'yes';
STUDY = std_selectdesign(STUDY, ALLEEG, STUDY.currentdesign);

com = sprintf('[STUDY ALLEEG] = pop_loadstudy(''filename'', ''%s'', ''filepath'', ''%s'');', STUDY.filename, STUDY.filepath);
