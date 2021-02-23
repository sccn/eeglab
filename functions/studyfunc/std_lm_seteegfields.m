% std_lm_seteegfields() - set limo fields in eeg sets
%
% Usage:
%   >>   EEG = std_lm_seteegfields(STUDY,index,'datatype','Channels','format', 'cell')
%
% Inputs:
%      STUDY    - studyset structure containing some or all files in EEG
%      index    - index of dataset in STUDY.datasetinfo
%
% Optional inputs:
%       format    - 'matrix' (default) or 'cell'
%       datatype  - 'channels' or 'components'
%
% Outputs:
%      EEG    - updated EEG structure
%
% See also:
%
% Authors: Cyril Pernet, The university of Edinburgh
%          Arnaud Delorme, SCCN
%          Ramon Martinez-Cancino, SCCN

% Copyright (C) 2015  Ramon Martinez-Cancino,INC, SCCN
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

function EEG = std_lm_seteegfields(STUDY,EEG,index,varargin)

%% Settings defaults
%  -----------------
try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            opt.(options{i}) = options{i+1};
        end
    else opt= []; end
catch
    error('limo_create_single_trials() error: calling convention {''key'', value, ... } error');
end

try, opt.format;           catch, opt.format        = 'cell';      end
try, opt.datatype;         catch, opt.datatype      = 'channels';  end

% Getting measureflags
% --------------------
measures   = {'erp','spec','ersp','timef','itc'};

% Getting prefix for channels/components
if strncmpi(opt.datatype,'chan',4)
    prefix = 'dat';
elseif strncmpi(opt.datatype,'comp',4)
    prefix = 'ica';
end

% Getting values
for i = 1:length(measures)
    opt.(measures{i}) = STUDY.etc.measureflags.([prefix measures{i}]);
end
% ---

fields = fieldnames(opt);
c = 1;
for i = 1: length(structfun(@numel,opt))
    if ~any ([isempty(eval(['opt.' fields{i}])),strcmp(fields(i),'format'), strcmp(fields(i),'datatype')])
    in_options{c}   = fields{i};
    in_options{c+1} = eval(['opt.' fields{i}]);
    c = c+2;
    end
end

%% Loading datase
%  --------------
path_tmp = rel2fullpath(STUDY.filepath,STUDY.datasetinfo(index).filepath); 
name = fullfile(path_tmp, STUDY.datasetinfo(index).subject);

%% Channels: update EEG.set file
%  -----------------------------
if strcmpi(opt.datatype,'channels')
    % DATERP 
    if strcmp(opt.erp,'on')
        if ~exist([name '.daterp'],'file')
            tmp = dir([name '*.daterp']);
            name = fullfile(tmp(1).folder,tmp(1).name);
            warning('couldn''t find a direct match between .set and .daterp\n loading %s',name)
            data = load('-mat',name);
        else
            data = load('-mat',[name '.daterp']);
        end
        EEG.etc.timeerp = data.times;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_daterp.mat'],'data'); clear data
            EEG.etc.datafiles.daterp = [name '_daterp.mat'];
            delete([name '.daterp']);
        else
            EEG.etc.datafiles.daterp = [name '.daterp'];
        end
    end
    
    % DATSPEC
    if strcmp(opt.spec,'on')
        if ~exist([name '.datspec'],'file')
            tmp = dir([name '*.datspec']);
            name = fullfile(tmp(1).folder,tmp(1).name);
            warning('couldn''t find a direct match between .set and .datspec\n loading %s',name)
            data = load('-mat',name);
        else
            data = load('-mat',[name '.datspec']);
        end
        EEG.etc.freqspec = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_datspec.mat'],'data'); clear data
            EEG.etc.datafiles.datspec = [name '_datspec.mat'];
            delete([name '.datspec']);
        else
            EEG.etc.datafiles.datspec = [name '.datspec'];
        end
    end
    % DATERSP    
    if strcmp(opt.timef,'on')
        if ~exist([name '.dattimef'],'file')
            tmp = dir([name '*.dattimef']);
            name = fullfile(tmp(1).folder,tmp(1).name);
            warning('couldn''t find a direct match between .set and .dattimef\n loading %s',name)
            data = load('-mat',name);
        else
            data = load('-mat',[name '.dattimef'],'times','freqs');
        end
        EEG.etc.timeersp = data.times;
        EEG.etc.freqersp = data.freqs;
        if strcmp(opt.format,'matrix')
            disp('reading single trials ersp, be patient ...')
            data = load('-mat',[name '.dattimef']);
            data = limo_struct2mat(data);
            save([name '_datersp.mat'],'data'); clear data
            EEG.etc.datafiles.datersp = [name '_datersp.mat'];
        else
            EEG.etc.datafiles.datersp = [name '.dattimef'];
        end
    end
    % DATITC
    if strcmp(opt.itc,'on')
        if ~exist([name '.datitc'],'file')
            tmp = dir([name '*.datitc']);
            name = fullfile(tmp(1).folder,tmp(1).name);
            warning('couldn''t find a direct match between .set and .datitc\n loading %s',name)
            data = load('-mat',name);
        else
            data = load('-mat',[name '.datitc']);
        end
        EEG.etc.timeitc = data.times;
        EEG.etc.freqitc = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_datitc.mat'],'data'); clear data
            EEG.etc.datafiles.datitc = [name '_datitc.mat'];
            delete([name '.datitc']);
        else
            EEG.etc.datafiles.datitc = [name '.datitc'];
        end
    end
end

%% Components: update EEG.set file
%  -------------------------------
if strcmpi(opt.datatype,'components')
    if strcmp(opt.erp,'on')
        if ~exist([name '.icaerp'],'file')
            tmp = dir([name '*.icaerp']);
            name = fullfile(tmp(1).folder,tmp(1).name);
            warning('couldn''t find a direct match between .set and .icaerp\n loading %s',name)
            data = load('-mat',name);
        else
            data = load('-mat',[name '.icaerp']);
        end
        EEG.etc.timeerp = data.times;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaerp.mat'],'data'); clear data
            EEG.etc.datafiles.daterp = [name '_icaerp.mat'];
            delete([name '.icaerp']);
        else
            EEG.etc.datafiles.icaerp = [name '.icaerp'];
        end
    end
    % ICAERP
    if strcmp(opt.spec,'on')
        if ~exist([name '.icaspec'],'file')
            tmp = dir([name '*.icaspec']);
            name = fullfile(tmp(1).folder,tmp(1).name);
            warning('couldn''t find a direct match between .set and .icaspec\n loading %s',name)
            data = load('-mat',name);
        else
            data = load('-mat',[name '.icaspec']);
        end
        EEG.etc.freqspec = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaspec.mat'],'data'); clear data
            EEG.etc.datafiles.datspec = [name '_icaspec.mat'];
            delete([name '.icaspec']);
        else
            EEG.etc.datafiles.icaspec = [name '.icaspec'];
        end
    end
    % ICAERSP    
    if strcmp(opt.timef,'on')
        if ~exist([name '.icatimef'],'file')
            tmp = dir([name '*.icatimef']);
            name = fullfile(tmp(1).folder,tmp(1).name);
            warning('couldn''t find a direct match between .set and .icatimef\n loading %s',name)
            data = load('-mat',name);
        else
            data = load('-mat',[name '.icatimef']);
        end
        EEG.etc.timeersp = data.times;
        EEG.etc.freqersp = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaersp.mat'],'data'); clear data
            EEG.etc.datafiles.datersp = [name '_icaersp.mat'];
            delete([name '.icatimef']);
        else
            EEG.etc.datafiles.icaersp = [name '.icatimef'];
        end
    end
    % ICAITC
    if strcmp(opt.itc,'on')
        if ~exist([name '.icaitc'],'file')
            tmp = dir([name '*.icaitc']);
            name = fullfile(tmp(1).folder,tmp(1).name);
            warning('couldn''t find a direct match between .set and ..icaitc\n loading %s',name)
            data = load('-mat',name);
        else
            data = load('-mat',[name '.icaitc']);
        end
        EEG.etc.timeitc = data.times;
        EEG.etc.freqitc = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaitc.mat'],'data'); clear data
            EEG.etc.datafiles.datitc = [name '_icaitc.mat'];
            delete([name '.icaitc']);
        else
            EEG.etc.datafiles.icaitc = [name '.icaitc'];
        end
    end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function file_fullpath = rel2fullpath(studypath,filepath)
% Return full path if 'filepath' is a relative path. The output format will
% fit the one of 'filepath'. That means that if 'filepath' is a cell array,
% then the output will a cell array too, and the same if is a string.

nit = 1; if iscell(filepath), nit = length(filepath);end

for i = 1:nit
    if iscell(filepath),pathtmp = filepath{i}; else pathtmp = filepath; end
    if strfind(pathtmp(end),filesep), pathtmp = pathtmp(1:end-1); end % Getting rid of filesep at the end
    if strfind(pathtmp(1:2),['.' filesep])
        if iscell(filepath),
            file_fullpath{i} = fullfile(studypath,pathtmp(3:end));
        else
            file_fullpath = fullfile(studypath,pathtmp(3:end));
        end
    else
        if iscell(filepath),
            file_fullpath{i} = pathtmp;
        else
            file_fullpath = pathtmp;
        end
    end
end
