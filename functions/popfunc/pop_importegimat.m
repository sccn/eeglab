% pop_importegimat() - import EGI Matlab segmented file
%
% Usage:    
%   >> EEG = pop_importegimat(filename, srate, latpoint0, getchanfile); 
%
% Inputs:
%   filename    - Matlab file name
%   srate       - sampling rate
%   latpoint0   - latency in sample ms of stimulus presentation.
%                 When data files are exported using Netstation, the user specify
%                 a time range (-100 ms to 500 ms for instance). In this 
%                 case, the latency of the stimulus is 100 (ms). Default is 0 (ms)
%   getchanfile - [0|1] automatically look up channel location file
%
% Output:
%   EEG        - EEGLAB dataset structure
%
% Authors: Arnaud Delorme, CERCO/UCSD, Jan 2010

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [EEG com] = pop_importegimat(filename, srate, latpoint0, lookupchanfile);

    if nargin < 3, latpoint0 = 0; end;
    if nargin < 4, lookupchanfile = 0; end;
    if nargin < 1 
        % ask user
        [filename, filepath] = uigetfile('*.mat', 'Choose a Matlab file from Netstation -- pop_importegimat()'); 
        if filename == 0 return; end;
        filename = fullfile(filepath, filename);

        promptstr    = { { 'style' 'text'       'string' 'Sampling rate (Hz)' } ...
                         { 'style' 'edit'       'string' '250' } ...
                         { 'style' 'text'       'string' 'Sample latency for stimulus (ms)' } ...
                         { 'style' 'edit'       'string' '0' } ...
                         { 'style' 'text'       'string' 'Look up channel location file (set=yes)' } ...
                         { 'style' 'checkbox'   'string' '' 'value' 0 } };
        geometry = { [2 1] [2 1] [2 1] };
        result       = inputgui( 'geometry', geometry, 'uilist', promptstr, ...
                                 'helpcom', 'pophelp(''pop_importegimat'')', ...
                                 'title', 'Import a Matlab file from Netstation -- pop_runica()');
                             
        if length(result) == 0 return; end;        
        srate = str2num(result{1});
        latpoint0 = str2num(result{2});
        if isempty(latpoint0), latpoint0 = 0; end;
        lookupchanfile = result{3};
    end;
   
    EEG = eeg_emptyset;
    fprintf('Reading EGI Matlab file %s\n', filename);
    tmpdata = load('-mat', filename);
    
    % get data types
    % --------------
    allfields = fieldnames(tmpdata);
    for index = 1:length(allfields)
        allfields{index} = allfields{index}(1:findstr(allfields{index}, 'Segment')-2);
    end;
    datatypes = unique(allfields);
    datatypes(cellfun(@isempty, datatypes)) = [];
    
    % read all data
    % -------------
    counttrial = 1;
    EEG.srate  = srate;
    latency = (latpoint0/1000)*EEG.srate+1;
    for index = 1:length(datatypes)
        tindex = 1;
        while isfield(tmpdata, sprintf('%s_Segment%d', datatypes{index}, tindex))
            datatrial = getfield(tmpdata, sprintf('%s_Segment%d', datatypes{index}, tindex));
            if counttrial == 1
                EEG.pnts = size(datatrial,2);
                EEG.data = repmat(single(0), [size(datatrial,1), size(datatrial,2), 1000]);
            end;
            EEG.data(:,:,counttrial) = datatrial;
            
            EEG.event(counttrial).type    = datatypes{index};
            EEG.event(counttrial).latency = latency;
            EEG.event(counttrial).epoch   = counttrial;
            
            counttrial = counttrial+1;
            tindex     = tindex + 1;
            latency = latency + EEG.pnts;
        end;
    end;
    fprintf('%d trials read\n', counttrial-1);
    EEG.data(:,:,counttrial:end) = [];
    
    EEG.setname  = filename(1:end-4);
    EEG.nbchan   = size(EEG.data,1);
    EEG.trials   = counttrial-1;
    if latpoint0 ~= 1
        EEG.xmin     = -latpoint0/1000;
    end;
    EEG = eeg_checkset(EEG);
    
    % channel location
    % ----------------
    if lookupchanfile
        EEG = readegilocs(EEG);
    end;
    
    com = sprintf('EEG = pop_importegimat(''%s'', %3.2f, %3.2f, %d);', filename, srate, latpoint0, lookupchanfile);
    