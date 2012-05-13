% pop_importpi() - EEGLAB plugin for importing Procom Infinity data files.
%
% Usage:
%   >> EEG = pop_importpi;
%
% Inputs:
%   filename   - [string] file name
%
% Optional inputs:
%   chanind    - [integer] channel indices (default all)
%
% Outputs:
%   EEG      - EEGLAB EEG strcuture
%
% Author: Arnaud Delorme, meditation research institute, Rishikesh, India

% Copyright (C) 2011 Arnaud Delorme, meditation research institute,
% Rishikesh, India
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

function [EEG com] = pop_importpi(filename, chanind_local)

    if nargin < 1
        [tmpf tmpp] = uigetfile('*.txt', 'Select Procom Infinity Text File');
        if tmpf(1) == 0, return; end;
        
        filename = fullfile(tmpp, tmpf);
        [chanlist srate] = getchanlist( filename);
        
        for index = 1:length(chanlist)
            chanlist2{index} = [ chanlist{index} ' at ' int2str(srate(index)) ' Hz' ];
        end;
        
        uilist = { { 'style' 'text' 'string' 'Import channel list' }
                   { 'style' 'listbox' 'string' chanlist2 'value' [1:length(chanlist)] 'max' 2} };
        res = inputgui('uilist', uilist, 'geometry', [1 1], 'geomvert', [1 3]);
        if isempty(res), return; end;
        chanind_local = res{1};
        
        if length(unique(srate(chanind_local))) > 1
            error('Cannot import data channels with different sampling rate');
        end;
        
    end;
        
    [chanlist srate] = getchanlist( filename);    
    if ~exist('chanind_local')
        chanind_local = 1:length(chanlist);
    end;
    
    % import data
    fid = fopen(filename, 'r');
    format = zeros(length(chanlist),2);
    format(1:length(chanlist),1) = '%';
    format(1:length(chanlist),2) = 'f';
    format = format';
    format = char(format(:)');
    %for index = 1:9, fgetl(fid); end;
    %realdata = fscanf(fid, '%f', [Inf]);
    realdata = textscan(fid, format, 'headerlines', 8, 'delimiter', ',');
    fclose(fid);
    
    % cop to EEG strcuture
    EEG = eeg_emptyset;
    for index = 1:length(chanind_local)
        EEG.data(index,:) = realdata{chanind_local(index)};
        EEG.chanlocs(index).labels = chanlist{chanind_local(index)};
    end;
    EEG.srate  = srate(1);
    EEG.nbchan = size(EEG.data,1);
    EEG = eeg_checkset(EEG);
    
    % output command
    if nargout == 2
        com = sprintf('EEG = pop_importpi(''%s'', [%s]);', filename, int2str(chanind_local));
    end;
 
function [chans srate ] = getchanlist(filename);

tmpa  = loadtxt(filename, 'skipline', 6, 'nlines', 2, 'delim', ', ', 'verbose', 'off');
chans = tmpa(1,1:3:end);
srate = [tmpa{2,:}];
