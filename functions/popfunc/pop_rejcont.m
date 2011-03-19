% pop_rejcont() - reject continuous portions of data based on spectrum
%                 thresholding. First, contiguous data epochs are extracted
%                 and a standard spectrum thresholding algorithm is
%                 applied. Regions of contiguous epochs larger than a
%                 specified size are then labeled as artifactual.
%
% Usage:
%   >> pop_rejcont( INEEG ) % pop-up interative window mode
%   >> [OUTEEG, selectedregions] = pop_rejcont( INEEG, 'key', 'val');
%
% Inputs:
%   INEEG      - input dataset
%
% Optional inputs:
%  'elecrange'   - [integer array] electrode indices {Default: all electrodes} 
%  'epochlength' - [float] epoch length in seconds {Default: 0.5 s}
%  'freqlimit'   - [min max] frequency range too consider for thresholding
%                  Default is [35 128] Hz.
%  'threshold'   - [float] frequency upper threshold in dB {Default: 10}
%  'contiguous'  - [integer] number of contiguous epochs necessary to 
%                  label a region as artifactual {Default: 4}
%  'addlength'   - [float] once a region of contiguous epochs has been labeled
%                  as artifact, additional trailing neighboring regions on
%                  each side may also be added {Default: 0.25}
%  'eegplot'     - ['on'|'off'] plot rejected portions of data in a eegplot
%                  window. Default is 'off'.
%
% Outputs:
%   OUTEEG          - output dataset with updated joint probability array
%   selectedregions - frames indices of rejected electrodes. Array of n x 2
%                     n being the number of regions and 2 for the beginning
%                     and end of each region.
%
% Author: Arnaud Delorme, CERCO, UPS/CNRS, 2009-
%
% See also: eegthresh()

% Copyright (C) 2009 Arnaud Delorme, CERCO, UPS/CNRS
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

function [EEG selectedregions com ] = pop_rejcont(EEG, varargin);

com = '';
if nargin < 1
    help pop_rejcont;
    return;
end;

if nargin < 2
    firstelec = 'EXG1';                % first non EEG channel
    % take all scalp electrodes
    % -------------------------
    if ~isempty(EEG.chanlocs)
        indelec = strmatch( firstelec, { EEG.chanlocs.labels });
        
        if isempty(indelec), elecrange = 1:EEG.nbchan;
        else                 elecrange = 1:(indelec-1);
        end;
    else
        elecrange = 1:EEG.nbchan;
    end;
    elecrange = deblank(vararg2str(elecrange));
    %elecrange = elecrange(2:end-1);
    
    promptstr = { 'Channel range' ...
                  'Frequency range (Hz)' ...
                  'Frequency threshold in dB' ...
                  'Epoch segment length (s)' ...
                  'Minimum number of contiguous epochs' ...
                  'Add trails before and after regions (s)' ...
                  };
    initstr = { elecrange '20 40' '10' '0.5' '4' '0.25' };
    result = inputdlg2(promptstr, 'Reject portions of continuous data - pop_rejcont', 1, initstr);
    if length( result ) == 0 return; end;
    
    options = { 'elecrange'     str2num(result{1}) ...
                'freqlimit'     str2double(result{2}) ...
                'threshold'     str2double(result{3}) ...
                'epochlength'   str2double(result{4}) ...
                'contiguous'    str2double(result{5}) ...
                'addlength'     str2double(result{6}) };
else
    options = varargin;
end;

opt = finputcheck(options, { 'threshold'   'real'   []    10;
                             'elecrange'   'real'   []    [1:EEG.nbchan];
                             'freqlimit'   'real'   []    [20 40];
                             'contiguous'  'real'   []    4;
                             'addlength'   'real'   []    0.25;
                             'eegplot'     'string' { 'on' 'off' } 'off';
                             'epochlength' 'real'   []    0.5 }, 'pop_rejcont');
     
%EEG.event = [];
grouplen  = opt.contiguous/2*opt.epochlength*EEG.srate+1; % maximum number of points for grouping regions
color     = [ 0 0.9 0];             % color of rejection window

% compute power spectrum
% ----------------------
NEWEEG = EEG;
% average reference 
% NEWEEG.data(opt.elecrange,:) = NEWEEG.data(opt.elecrange,:)-repmat(mean(NEWEEG.data(opt.elecrange,:),1), [length(opt.elecrange) 1]);

[TMPNEWEEG]= eeg_regepochs(NEWEEG, 0.25, [0 opt.epochlength], NaN);
%[TMPNEWEEG indices] = pop_rejspec(TMPNEWEEG, 1, [1:64], -100, 15, 30, 45, 0, 0);
%rejepoch = find(indices);
tmp   = fft(TMPNEWEEG.data, [], 2);
freqs = linspace(0, TMPNEWEEG.srate/2, size(tmp,2)/2);
freqs = freqs(2:end); % remove DC (match the output of PSD)
tmp   = tmp(:,2:size(tmp,2)/2,:);
warning('off', 'MATLAB:log:logOfZero');
tmpspec  = 10*log10(abs(tmp).^2);  
warning('on', 'MATLAB:log:logOfZero');
tmpspec  = tmpspec - repmat( mean(tmpspec,3), [1 1 TMPNEWEEG.trials]);
specdata = tmpspec;

% apply threshold to average of all electrodes
% --------------------------------------------
%[I1 Irej NS Erej] = eegthresh( mean(specdata(opt.elecrange, :, :), 1), size(specdata,2), 1:length(opt.elecrange), -100, 15, [freqs(1) freqs(end)], 30, 45);
if length(opt.threshold) == 1, opt.threshold = [ -100 opt.threshold ]; end;
[I1 rejepoch NS Erej] = eegthresh( mean(specdata(opt.elecrange, :, :), 1), size(specdata,2), 1, opt.threshold(1), opt.threshold(2), [freqs(1) freqs(end)], opt.freqlimit(1), opt.freqlimit(2));
fprintf('%d regions selected for rejection\n', length(rejepoch));

% build the winrej array for eegplot
% ----------------------------------
winrej   = [];
if ~isempty(find(cellfun(@isempty, { TMPNEWEEG.event.epoch }) == 1))
    error('Some events are not associated with any epoch');
end;
allepoch = [ TMPNEWEEG.event.epoch ];
for index = 1:length(rejepoch)
    eventepoch = find( rejepoch(index) == allepoch );
    if strcmpi(TMPNEWEEG.event(eventepoch(1)).type, 'X')
        urevent = TMPNEWEEG.event(eventepoch(1)).urevent;
        lat     = TMPNEWEEG.urevent(urevent).latency;
        %if length(
        winrej  = [ winrej; lat lat+opt.epochlength*NEWEEG.srate-1 color ]; %Erej(:,index)'];
    else
        error('Wrong type for epoch');
    end;
end;
winrej(:,6:6+length(opt.elecrange)-1) = 0;

% remove isolated regions and merge others
% ----------------------------------------
merged = 0;
isolated = 0;
for index = size(winrej,1):-1:1
    if size(winrej,1) >= index && winrej(index,2) - winrej(index,1) > grouplen, winrej(index,:) = []; isolated = isolated + 1;
    elseif index == 1 && size(winrej,1) > 1 && winrej(index+1,1) - winrej(index,2) > grouplen, winrej(index,:) = []; isolated = isolated + 1;
    elseif index == size(winrej,1) && size(winrej,1) > 1 && winrej(index,1) - winrej(index-1,2) > grouplen, winrej(index,:) = []; isolated = isolated + 1;
    elseif index > 1 && size(winrej,1) > 1 && index < size(winrej,1) && winrej(index+1,1) - winrej(index,2) > grouplen && ...
            winrej(index,1) - winrej(index-1,2) > grouplen
        winrej(index,:) = [];
        isolated = isolated + 1;
    elseif index < size(winrej,1) && size(winrej,1) > 1 && winrej(index+1,1) - winrej(index,2) <= grouplen
        winrej(index,2) = winrej(index+1,2);
        winrej(index+1,:) = [];
        merged = merged + 1;
    end;
end;
fprintf('%d regions merged\n', merged);
fprintf('%d regions removed\n', isolated);

% add time before and after each region
% -------------------------------------
for index = 1:size(winrej,1)
    winrej(index,1) = max(1,         winrej(index,1)-opt.addlength);
    winrej(index,2) = min(EEG.pnts,  winrej(index,2)+opt.addlength);
end;

% plot result
% -----------
selectedregions = winrej(:,1:2);
command = 'EEG = pop_select(EEG, ''nopoint'', TMPREJ(:,1:2)); [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, ''study'', ~isempty(STUDY)+0); eeglab redraw';
if nargin < 2 || strcmpi(opt.eegplot, 'on')
    eegplot(NEWEEG.data(opt.elecrange,:), 'srate', NEWEEG.srate, 'winrej', winrej, 'command', command, 'events', EEG.event);
else
    NEWEEG = pop_select(EEG, 'nopoint', round(selectedregions));
end;
EEG = NEWEEG;
com = sprintf('EEG = pop_rejcont(EEG, %s);', vararg2str(options));
