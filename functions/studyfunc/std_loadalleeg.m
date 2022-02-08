% std_loadalleeg() - constructs an ALLEEG structure, given the paths and file names 
%                    of all the EEG datasets that will be loaded into the ALLEEG 
%                    structure. The EEG datasets may be loaded without their EEG.data 
%                    (see the pop_editoptions() function), so many datasets can be 
%                    loaded simultaneously. The loaded EEG datasets have dataset 
%                    information and a (filename) pointer to the data. 
% Usage:    
%                  % Load sseveral EEG datasets into an ALLEEG structure.
%                  >> ALLEEG = std_loadalleeg(paths,datasets) ;  
%                  >> ALLEEG = std_loadalleeg(STUDY) ;  
% Inputs:
%   paths      - [cell array of strings] cell array with all the datasets paths. 
%   datasets   - [cell array of strings] cell array with all the datasets file names. 
%
% Output:
%   ALLEEG     - an EEGLAB data structure, which holds the loaded EEG sets  
%                 (can also be one EEG set).
%  Example: 
%          >> paths = {'/home/eeglab/data/sub1/','/home/eeglab/data/sub2/', ...  
%          >>          '/home/eeglab/data/sub3/', '/home/eeglab/data/sub6/'};
%          >> datasets = { 'visattS1', 'visattS2', 'visattS3', 'visattS4'};
%          >> ALLEEG = std_loadalleeg(paths,datasets) ; 
%                
% See also: pop_loadstudy(), pop_study()
%
% Authors: Hilit Serby, Arnaud Delorme, SCCN, INC, UCSD, October , 2004

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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

function ALLEEG = std_loadalleeg(varargin)

    if nargin < 1
        help std_loadalleeg;
        return;
    end

    genpath = '';
    oldgenpath = '';
    if isstruct(varargin{1})
        datasets = {varargin{1}.datasetinfo.filename};
        try,
            paths = {varargin{1}.datasetinfo.filepath};
        catch,
            paths = cell(1,length(datasets));
            paths(:) = { '' };
        end
        genpath = varargin{1}.filepath;
        if isfield(varargin{1}.etc, 'oldfilepath')
            oldgenpath = varargin{1}.etc.oldfilepath;
        end
    else
        paths = varargin{1};
        if nargin > 1
            datasets = varargin{2};
        else
            datasets = paths;
            paths    = cell(size(datasets));
        end
    end
   
    set = 1;
    ALLEEG = [];
    
    eeglab_options;
    
    % read datasets
    % -------------
    comp = computer;
    warnfold = 'off';
    if length(oldgenpath) > 1 && oldgenpath(2) == ':' && ~strcmpi(comp(1:2), 'PC')
        oldgenpath = [ filesep oldgenpath(4:end) ];
        oldgenpath(find(oldgenpath == '\')) = filesep;
    end
    
    for dset = 1:length(paths)
        if ~isempty(paths{dset})
            if paths{dset}(2) == ':' && ~strcmpi(comp(1:2), 'PC') 
                paths{dset} = [ filesep paths{dset}(4:end) ];
                paths{dset}(find(paths{dset} == '\')) = filesep;
            end
        end
        
        [sub2 sub1] = fileparts(char(paths{dset}));
        [sub3 sub2] = fileparts(sub2);
        
        % priority is given to relative path of the STUDY if STUDY has moved
        if ~isempty(oldgenpath) && oldgenpath(end) == filesep
             oldgenpath(end) = [];
        end
        if ~isequal(genpath, oldgenpath) && ~isempty(oldgenpath)
            disp('Warning: STUDY moved since last saved, trying to load data files using relative path');
            if  ~isempty(strfind(char(paths{dset}), oldgenpath))
                relativePath = char(paths{dset}(length(oldgenpath)+1:end));
                relativePath = fullfile(genpath, relativePath);
            else
                disp('Important warning: relative path cannot calculated, make sure the correct data files are loaded');
                relativePath = char(paths{dset});
            end;   
            
            % fix issue when datasets are in a parent folder of the STUDY
            if dset == 1
                indCommon = 1;
                while indCommon <= length(oldgenpath) && indCommon <= length(paths{1}) && paths{1}(indCommon) == oldgenpath(indCommon)
                    indCommon = indCommon+1;
                end
                indCommon = indCommon-1;
                if indCommon > 1 && indCommon < length(genpath) % do not change path if nothing in common between the two paths
                    genpath(indCommon-length(oldgenpath)+length(genpath)+1:end) = [];
                end
            end
        else
            relativePath = char(paths{dset});
        end
        
        % load data files
        if exist(fullfile(relativePath, datasets{dset})) == 2
            EEG = pop_loadset('filename', datasets{dset}, 'filepath', relativePath, 'loadmode', 'info', 'check', 'off');
        elseif exist(fullfile(char(paths{dset}), datasets{dset})) == 2
            EEG = pop_loadset('filename', datasets{dset}, 'filepath', char(paths{dset}), 'loadmode', 'info', 'check', 'off');
        elseif exist( fullfile(genpath, datasets{dset})) == 2    
            [tmpp tmpf ext] = fileparts(fullfile(genpath, datasets{dset}));
            EEG = pop_loadset('filename', [tmpf ext], 'filepath',tmpp, 'loadmode', 'info', 'check', 'off');
            warnfold = 'on';
        elseif exist( fullfile(genpath, sub1, datasets{dset})) == 2    
            [tmpp tmpf ext] = fileparts(fullfile(genpath, sub1, datasets{dset}));
            EEG = pop_loadset('filename', [tmpf ext], 'filepath',tmpp, 'loadmode', 'info', 'check', 'off');
            warnfold = 'on';
        elseif exist( fullfile(genpath, sub2, datasets{dset}))  == 2   
            [tmpp tmpf ext] = fileparts(fullfile(genpath, sub2, datasets{dset}));
            EEG = pop_loadset('filename', [tmpf ext], 'filepath',tmpp, 'loadmode', 'info', 'check', 'off');
            warnfold = 'on';
        elseif exist( fullfile(genpath, sub2, sub1, datasets{dset}))  == 2   
            [tmpp tmpf ext] = fileparts(fullfile(genpath, sub2, sub1, datasets{dset}));
            EEG = pop_loadset('filename', [tmpf ext], 'filepath',tmpp, 'loadmode', 'info', 'check', 'off');
            warnfold = 'on';
        elseif exist(lower(fullfile(char(paths{dset}), datasets{dset}))) == 2   
            EEG = pop_loadset('filename', lower(datasets{dset}), 'filepath',lower(char(paths{dset})), 'loadmode', 'info', 'check', 'off');
        elseif exist(fullfile(pwd, datasets{dset})) == 2
            EEG = pop_loadset('filename', lower(datasets{dset}), 'filepath',pwd, 'loadmode', 'info', 'check', 'off');
        else
            txt = [ sprintf('The dataset %s is missing\n', datasets{dset}) 10 ...
                          'Is it possible that it might have been deleted?' 10 ...
                          'If this is the case, re-create the STUDY using the remaining datasets' ];
            error(txt);
        end
        
        EEG = eeg_checkset(EEG);
        if ~option_storedisk
            EEG = eeg_checkset(EEG, 'loaddata');
        elseif ~ischar(EEG.data)
            EEG.data   = 'in set file';
            EEG.icaact = [];
        end
        
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 0, 'notext');
    end
    ALLEEG = eeg_checkset(ALLEEG);

    if strcmpi(warnfold, 'on') && ~strcmpi(pwd, genpath) && ~isempty(genpath)
        disp('Changing current path to STUDY path...');
        cd(genpath);
    end
    if strcmpi(warnfold, 'on') 
        disp('This STUDY has a relative path set for the datasets');
        disp('so the current path MUST remain the STUDY path');
    end
        
