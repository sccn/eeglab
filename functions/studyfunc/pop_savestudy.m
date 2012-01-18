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
%
% Note: the parameter EEG is currenlty not being used. In the future, this function
%       will check if any of the datasets of the study have been modified and
%       have to be resaved.
%
% Authors: Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, September 2005

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

function [STUDY, EEG, com] = pop_savestudy(STUDY, EEG, varargin);

com = '';
if nargin < 1
	help pop_savestudy;
	return;
end;
if isempty(STUDY)  , error('pop_savestudy(): cannot save empty STUDY'); end;
if length(STUDY) >1, error('pop_savestudy(): cannot save multiple STUDY sets'); end;

% backward compatibility
% ----------------------
if nargin > 1 
    if isstr(EEG)
        options = { EEG varargin{:} };
    else
        options = varargin;
    end;
end;

if nargin < 3
    % pop up window to ask for file type
    % ----------------------------------
    [filename, filepath] = uiputfile2('*.study', ...
                    'Save STUDY with .study extension -- pop_savestudy()'); 
    if isequal(filename,0), return; end;
    if ~strncmp(filename(end-5:end), '.study',6)
        if isempty(strfind(filename,'.'))
            filename = [filename '.study'];
        else
            filename = [filename(1:strfind(filename,'.')-1) '.study'];
        end
    end
    options = { 'filename' filename 'filepath' filepath };
end

% decoding parameters
% -------------------
g = finputcheck(options,  { 'filename'   'string'   []     STUDY.filename;
                            'filepath'   'string'   []     STUDY.filepath;
                            'savemode'   'string'   { 'standard','resave' } 'standard' });
if isstr(g), error(g); end;

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
    end;
    if isfield(STUDY.changrp, fields{fInd})
        STUDY.cluster = rmfield(STUDY.cluster, fields{fInd});
    end;
end;    

% resave mode
% -----------
STUDY.saved = 'yes';
if strcmpi(g.savemode, 'resave')
    disp('Re-saving study file');
    g.filename = STUDY.filename;
    g.filepath = STUDY.filepath;
end;

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
v = version;
if v(1) > '6'
    save('-v6' , STUDYfile, 'STUDY');
else
    save('-mat', STUDYfile, 'STUDY');
end;
STUDY = STUDYTMP;

% history
% -------
com = sprintf('[STUDY EEG] = pop_savestudy( %s, %s, %s);', inputname(1), inputname(2), vararg2str(options));
