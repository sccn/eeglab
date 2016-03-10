% std_lm_seteegfields() - set limo fields in eeg sets
%
% Usage:
%   >>   EEG = std_lm_seteegfields(STUDY,index,'datatype','Channels','format', 'cell')
%
% Inputs:
%      STUDY    - studyset structure containing some or all files in ALLEEG
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
% Author: Cyril Pernet, The university of Edinburgh
%         Arnaud Delorme, SCCN
%         Ramon Martinez-Cancino, SCCN
%
% Copyright (C) 2015  Ramon Martinez-Cancino,INC, SCCN
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

function EEG = std_lm_seteegfields(STUDY,index,varargin)

%% Settings defaults
%  -----------------
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            opt.(options{i}) = options{i+1};
        end
    else opt= []; end;
catch
    error('limo_create_single_trials() error: calling convention {''key'', value, ... } error');
end;

try, opt.format;           catch, opt.format        = 'cell';      end;
try, opt.datatype;         catch, opt.datatype      = 'channels';  end;

% Getting measureflags
% --------------------
measures   = {'erp','spec','ersp','itc'};

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
EEG = fullfile(path_tmp,STUDY.datasetinfo(index).filename);
ALLEEG = pop_loadset(EEG);
name = fullfile(path_tmp, STUDY.datasetinfo(index).subject);

%% Channels: update EEG.set file
%  -----------------------------
if strcmpi(opt.datatype,'channels')
    % DATERP 
    if strcmp(opt.erp,'on')
        data = load('-mat',[name '.daterp']);
        ALLEEG.etc.timeerp = data.times;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_daterp.mat'],'data'); clear data
            ALLEEG.etc.datafiles.daterp = [name '_daterp.mat'];
            delete([name '.daterp']);
        else
            ALLEEG.etc.datafiles.daterp = [name '.daterp'];
        end
    end
    % DATSPEC
    if strcmp(opt.spec,'on')
        data = load('-mat',[name '.datspec']);
        ALLEEG.etc.freqspec = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_datspec.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datspec = [name '_datspec.mat'];
            delete([name '.datspec']);
        else
            ALLEEG.etc.datafiles.datspec = [name '.datspec'];
        end
    end
    % DATERSP    
    if strcmp(opt.ersp,'on')
        disp('reading single trials ersp, be patient ..')
        data = load('-mat',[name '.dattimef']);
        ALLEEG.etc.timeersp = data.times;
        ALLEEG.etc.freqersp = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_datersp.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datersp = [name '_datersp.mat'];
            delete([name '.dattimef']);
        else
            ALLEEG.etc.datafiles.datersp = [name '.dattimef'];
        end
    end
    % DATITC
    if strcmp(opt.itc,'on')
        data = load('-mat',[name '.datitc']);
        ALLEEG.etc.timeitc = data.times;
        ALLEEG.etc.freqitc = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_datitc.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datitc = [name '_datitc.mat'];
            delete([name '.datitc']);
        else
            ALLEEG.etc.datafiles.datitc = [name '.datitc'];
        end
    end
end

%% Components: update EEG.set file
%  -------------------------------
if strcmpi(opt.datatype,'components')
    if strcmp(opt.erp,'on')
        data = load('-mat',[name '.icaerp']);
        ALLEEG.etc.timeerp = data.times;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaerp.mat'],'data'); clear data
            ALLEEG.etc.datafiles.daterp = [name '_icaerp.mat'];
            delete([name '.icaerp']);
        else
            ALLEEG.etc.datafiles.icaerp = [name '.icaerp'];
        end
    end
    % ICAERP
    if strcmp(opt.spec,'on')
        data = load('-mat',[name '.icaspec']);
        ALLEEG.etc.freqspec = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaspec.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datspec = [name '_icaspec.mat'];
            delete([name '.icaspec']);
        else
            ALLEEG.etc.datafiles.icaspec = [name '.icaspec'];
        end
    end
    % ICAERSP    
    if strcmp(opt.ersp,'on')
        data = load('-mat',[name '.icatimef']);
        ALLEEG.etc.timeersp = data.times;
        ALLEEG.etc.freqersp = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaersp.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datersp = [name '_icaersp.mat'];
            delete([name '.icatimef']);
        else
            ALLEEG.etc.datafiles.icaersp = [name '.icatimef'];
        end
    end
    % ICAITC
    if strcmp(opt.itc,'on')
        data = load('-mat',[name '.icaitc']);
        ALLEEG.etc.timeitc = data.times;
        ALLEEG.etc.freqitc = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaitc.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datitc = [name '_icaitc.mat'];
            delete([name '.icaitc']);
        else
            ALLEEG.etc.datafiles.icaitc = [name '.icaitc'];
        end
    end
end
%%
EEG = ALLEEG; clear ALLEEG;
pop_saveset(EEG, 'filename', EEG.filename, 'filepath',EEG.filepath,'savemode' ,'twofiles');
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
    if strfind(pathtmp,['.' filesep])
        if iscell(filepath),
            file_fullpath{i} = fullfile(studypath,pathtmp(3:end));
        else
            file_fullpath = fullfile(studypath,pathtmp(3:end));
        end
    else
        file_fullpath = pathtmp;
    end
end
end
