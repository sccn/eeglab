% load_alleeg() - This function constructs an ALLEEG structure.
% The function takes the paths and file names of all the EEG datasets that
% will be loaded into the ALLEEG structure. 
% The EEG datasets are loaded without their data and icaact fields to save
% memory space, so many datasets can be loaded simultaneously.
% The loaded EEG datasets have the set information and a pointer to the
% data. The datasets are needed to be saved before hand in this format
% using pop_saveset with the input argument 'savemode' set to 'twofiles' (see
% pop_saveset for details).
%
% Usage:    
%   >> ALLEEG = load_alleeg(paths,datasets) ;  
%   >> ALLEEG = load_alleeg(STUDY) ;  
%   The function loads several EEG datasets into an ALLEEG structure.
%
%
% Inputs:
%   paths           - [cell array of strings] cell array with all the datasets paths. 
%   datasets      - [cell array of strings] cell array with all the datasets file names. 
%
% Output:
%   ALLEEG     - an EEGLAB data structure, which holds the loaded EEG sets  
%                      (can also be one EEG set). The loaded EEG datasets
%                      don't include the data and ICA activation (EEG.data, EEG.icaact), 
%                      instead EEG.data holds the pointer to the floating file that 
%                      holds the dataset EEG data. The datasets are needed to be saved 
%                      before hand in this format using pop_saveset with the input 
%                      argument 'savemode' set to 'twofiles' (see pop_saveset for details).
%
%  example: 
%      paths = {'/home/eeglab/data/sub1/','/home/eeglab/data/sub2/', ...  
%                     '/home/eeglab/data/sub3/', '/home/eeglab/data/sub6/'};
%      datasets = { 'visattS1', 'visattS2', 'visattS3', 'visattS4'};
%      ALLEEG = load_alleeg(paths,datasets) ; 
%                
%  See also  pop_loadstudy, pop_createstudy, create_study          
%
% Authors:  Hilit Serby, Arnaud Delorme, SCCN, INC, UCSD, October , 2004

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

function ALLEEG = load_alleeg(varargin)

    if nargin < 1
        help load_alleeg;
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
    
    eeg_optionsbackup;
    eeg_options;
    
    % read datasets
    % -------------
    for dset = 1:length(paths)
        if exist(fullfile(paths{dset}, datasets{dset}))    
            EEG = pop_loadset(datasets{dset}, paths{dset}, 'info');
        elseif exist( fullfile(genpath, paths{dset}, datasets{dset}))    
            [tmpp tmpf ext] = fileparts(fullfile(genpath, paths{dset}, datasets{dset}));
            EEG = pop_loadset([tmpf ext], tmpp, 'info');
        else
            error(sprintf('Dataset ''%s'' not found', fullfile(paths{dset}, datasets{dset})));
        end;
        
        if ~option_storedisk
            EEG = eeg_checkset(EEG, 'loaddata');
        elseif ~isstr(EEG.data)
            EEG.data   = 'in set file';
            EEG.icaact = [];
        end;
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 0);
    end
