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

function ALLEEG = std_loadalleeg(varargin)

    if nargin < 1
        help std_loadalleeg;
        return;
    end;

    if isstruct(varargin{1})
        datasets = {varargin{1}.datasetinfo.filename};
        try,
            paths = {varargin{1}.datasetinfo.filepath};
        catch,
            paths = cell(1,length(datasets));
            paths(:) = { '' };
        end;
        genpath = varargin{1}.filepath;
    else
        genpath = '';
        paths = varargin{1};
        if nargin > 1
            datasets = varargin{2};
        else
            datasets = paths;
            paths    = cell(size(datasets));
        end;
    end
    
    set = 1;
    ALLEEG = [];
    
    eeglab_options;
    
    % read datasets
    % -------------
    warnfold = 'off';
    for dset = 1:length(paths)
        if ~isempty(paths{dset})
            comp = computer;
            if paths{dset}(2) == ':' & ~strcmpi(comp(1:2), 'PC') 
                paths{dset} = [ filesep paths{dset}(4:end) ];
                paths{dset}(find(paths{dset} == '\')) = filesep;
            end;
        end;
        
        [sub2 sub1] = fileparts(char(paths{dset}));
        [sub3 sub2] = fileparts(sub2);
        
        if exist(fullfile(char(paths{dset}), datasets{dset})) == 2
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
        else
            txt = [ sprintf('The dataset %s is missing\n', datasets{dset}) 10 ...
                          'Is it possible that it might have been deleted?' 10 ...
                          'If this is the case, re-create the STUDY using the remaining datasets' ];
            error(txt);
        end;
        
        EEG = eeg_checkset(EEG);
        if ~option_storedisk
            EEG = eeg_checkset(EEG, 'loaddata');
        elseif ~isstr(EEG.data)
            EEG.data   = 'in set file';
            EEG.icaact = [];
        end;
        
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 0, 'notext');
    end
    ALLEEG = eeg_checkset(ALLEEG);

    if strcmpi(warnfold, 'on') & ~strcmpi(pwd, genpath)
        disp('Changing current path to STUDY path...');
        cd(genpath);
    end;
    if strcmpi(warnfold, 'on') 
        disp('This STUDY has a relative path set for the datasets');
        disp('so the current path MUST remain the STUDY path');
    end;
        
