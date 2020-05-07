% pop_savestudy() - save a STUDY structure to a disk file
%
% Usage:
%   >> STUDY = pop_savestudy( STUDY, EEG ); % pop up and interactive window 
%   >> STUDY = pop_savestudy( STUDY, EEG, 'key', 'val', ...); % no pop-up
%                                              
% Inputs:
%   STUDY      - STUDY structure. 
%   EEG        - Array of datasets contained in the study.
%
% Optional inputs:
%   'filename' - [string] name of the STUDY file {default: STUDY.filename}
%   'filepath' - [string] path of the STUDY file {default: STUDY.filepath}
%   'savemode' - ['resave'|'standard'] in resave mode, the file name in
%                the study is being used to resave it.
%   'resavedatasets' - ['on'|'off'] resave STUDY datasets if necessary.
%                Default is 'off'.
%
% Note: the parameter EEG is currenlty not being used. In the future, this function
%       will check if any of the datasets of the study have been modified and
%       have to be resaved.
%
% Authors: Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, September 2005

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

function [STUDY, EEG, com] = pop_savestudy(STUDY, EEG, varargin);

com = '';
if nargin < 1
	help pop_savestudy;
	return;
end
if isempty(STUDY)  , error('pop_savestudy(): cannot save empty STUDY'); end
if length(STUDY) >1, error('pop_savestudy(): cannot save multiple STUDY sets'); end

% backward compatibility
% ----------------------
if nargin > 1 
    if ischar(EEG)
        options = { EEG varargin{:} };
    else
        options = varargin;
    end
end

if nargin < 3

    % pop up window to ask for file type
    % ----------------------------------
    [filename, filepath] = uiputfile2('*.study', ...
                    'Save STUDY with .study extension -- pop_savestudy()'); 
    if isequal(filename,0), return; end
    if ~strncmp(filename(end-5:end), '.study',6)
        if isempty(strfind(filename,'.'))
            filename = [filename '.study'];
        else
            filename = [filename(1:strfind(filename,'.')-1) '.study'];
        end
    end
        
    options = { options{:} 'filename' filename 'filepath' filepath };
end

% check if datasets have been modified and need to be saved
% ---------------------------------------------------------
if ~isempty(EEG)
    EEG = eeg_checkset(EEG);
    if any(cellfun(@(x)isequal(x, 'no'), { EEG.saved }))
        if nargin < 3 || (length(options) == 2 && strcmpi(options{2}, 'resavegui'))
            res = questdlg2('Some datasets have been modified, do you want to resave them?', ...
                'Saving datasets', ...
                'No', 'Yes', 'Yes');
            if strcmpi(res, 'yes')
                options = { options{:} 'resavedatasets' 'on' };
            end
        else
            disp('Some datasets have been modified, use option ''resavedatasets'' to resave them');
        end
    end
end

% decoding parameters
% -------------------
g = finputcheck(options,  { 'filename'   'string'   []     STUDY.filename;
                            'filepath'   'string'   []     STUDY.filepath;
                            'resavedatasets'   'string'   { 'on' 'off' } 'off';
                            'savemode'   'string'   { 'standard','resave','resavegui' } 'standard' });
if ischar(g), error(g); end
if isempty(STUDY.filename) && isempty(g.filename)
    error('File name required to save the study');
end

% resave datasets
% ---------------
if strcmpi(g.resavedatasets, 'on')
    EEG = pop_saveset(EEG, 'savemode', 'resave');
end

% fields to remove
% ----------------
fields = { 'erptimes'  'erpdata' ...
           'specfreqs' 'specdata'  ...
           'erspdata'  'ersptimes' 'erspfreqs' 'erspdatatrials' 'erspsubjinds' 'erspbase'    'ersptrialinfo' ...
           'itcdata'   'itcfreqs' 'itctimes' ...
           'erpimdata' 'erpimevents' 'erpimtrials' 'erpimtimes' };
for fInd = 1:length(fields)
    if isfield(STUDY.changrp, fields{fInd})
        STUDY.changrp = rmfield(STUDY.changrp, fields{fInd});
    end
    if isfield(STUDY.changrp, fields{fInd})
        STUDY.cluster = rmfield(STUDY.cluster, fields{fInd});
    end
end;    

% resave mode
% -----------
STUDY.saved = 'yes';
if strcmpi(g.savemode, 'resave')
    disp('Re-saving study file');
    g.filename = STUDY.filename;
    g.filepath = STUDY.filepath;
end

if isempty(g.filename)
    disp('pop_savestudy(): no STUDY filename: make sure the STUDY has a filename');
    return;
end
if ~strncmp(g.filename(end-5:end), '.study',6)
    if isempty(strfind(g.filename,'.'))
        g.filename = [g.filename '.study'];
    else
        g.filename = [g.filename(1:strfind(g.filename,'.')-1) '.study'];
    end
end

%    [filepath filenamenoext ext] = fileparts(varargin{1});
%    filename = [filenamenoext '.study']; % make sure a .study extension
    
STUDY.filepath = g.filepath;
STUDY.filename = g.filename;
STUDYfile = fullfile(STUDY.filepath,STUDY.filename);
STUDYTMP = STUDY;
STUDY = std_rmalldatafields(STUDY);
eeglab_options;
if option_saveversion6, save('-v6' , STUDYfile, 'STUDY');
else                    save('-v7.3' , STUDYfile, 'STUDY');
end
STUDY = STUDYTMP;

% history
% -------
com = sprintf('[STUDY EEG] = pop_savestudy( STUDY, EEG, %s);', vararg2str(options));
