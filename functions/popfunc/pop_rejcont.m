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
%  'overlap'     - [float] epoch overlap in seconds {Default: 0.25 s}
%  'freqlimit'   - [min max] frequency range too consider for thresholding
%                  Default is [35 128] Hz.
%  'mode'        - ['max'|'mean'] average power or take the max in the 
%                  selected frequency range. Default is 'max'.
%  'correct'     - ['remove'|'blank'] type of correction. Default is to
%                  'remove' the bad portion. 'Blank' put to 0 the selected
%                   electrodes.
%  'threshold'   - [float] frequency upper threshold in dB {Default: 10}
%  'contiguous'  - [integer] number of contiguous epochs necessary to 
%                  label a region as artifactual {Default: 4 }
%  'addlength'   - [float] once a region of contiguous epochs has been labeled
%                  as artifact, additional trailing neighboring regions on
%                  each side may also be added {Default: 0.25 s}
%  'eegplot'     - ['on'|'off'] plot rejected portions of data in a eegplot
%                  window. Default is 'off'.
%  'onlyreturnselection'  - ['on'|'off'] this option when set to 'on' only
%                  return the selected regions and does not remove them 
%                  from the datasets. This allow to perform quick
%                  optimization of the rejected portions of data.
%  'precompstruct' - [struct] structure containing precomputed spectrum (see
%                  Outputs) to be used instead of computing the spectrum.
%  'verbose'       - ['on'|'off'] display information. Default is 'off'.
%  'taper'         - ['none'|'hamming'] taper to use before FFT. Default is
%                    'none' for backward compatibility but 'hamming' is
%                    recommended.
%
% Outputs:
%   OUTEEG          - output dataset with region removed
%   selectedregions - frames indices of rejected electrodes. Array of n x 2
%                     n being the number of regions and 2 for the beginning
%                     and end of each region.
%   precompstruct   - structure containing precomputed data. This structure
%                     contains the spectrum, the frequencies and the EEGLAB
%                     dataset used as input with epochs extracted.
%
% Author: Arnaud Delorme, CERCO, UPS/CNRS, 2009-
%
% Example:
% EEG = pop_rejcont(EEG, 'elecrange',[1:32] ,'freqlimit',[20 40] ,'threshold',...
%    10,'epochlength',0.5,'contiguous',4,'addlength',0.25, 'taper', 'hamming');
% 
% See also: eegthresh()

% Copyright (C) 2009 Arnaud Delorme, CERCO, UPS/CNRS
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

function [EEG, selectedregions, precompstruct, com ] = pop_rejcont(EEG, varargin);

com = '';
selectedregions = [];
precompstruct = [];
if nargin < 1
    help pop_rejcont;
    return;
end

if nargin < 2
    firstelec = 'EXG1';                % first non EEG channel
    % take all scalp electrodes
    % -------------------------
    if ~isempty(EEG.chanlocs)
        tmpchanlocs = EEG.chanlocs;
        indelec = strmatch( firstelec, { tmpchanlocs.labels });
        
        if isempty(indelec), elecrange = 1:EEG.nbchan;
        else                 elecrange = 1:(indelec-1);
        end
    else
        elecrange = 1:EEG.nbchan;
    end
    elecrange = deblank(vararg2str(elecrange));
    %elecrange = elecrange(2:end-1);
    
%     promptstr = { 'Channel range' ...
%                   'Frequency range (Hz)' ...
%                   'Frequency threshold in dB' ...
%                   'Epoch segment length (s)' ...
%                   'Minimum number of contiguous epochs' ...
%                   'Add trails before and after regions (s)' ...
%                   };
%     initstr = { elecrange '20 40' '10' '0.5' '4' '0.25' };
%     result = inputdlg2(promptstr, 'Reject portions of continuous data - pop_rejcont', 1, initstr);
    uilist = { { 'style' 'text' 'string' 'Channel range' } ...
               { 'style' 'edit' 'string' elecrange } ...
               { 'style' 'text' 'string' 'Frequency range (Hz)' } ...
               { 'style' 'edit' 'string' '20 40' } ...
               { 'style' 'text' 'string' 'Frequency threshold in dB' } ...
               { 'style' 'edit' 'string' '10' } ...
               { 'style' 'text' 'string' 'Epoch segment length (s)' } ...
               { 'style' 'edit' 'string' '0.5' } ...
               { 'style' 'text' 'string' 'Minimum number of contiguous epochs' } ...
               { 'style' 'edit' 'string' '4' } ...
               { 'style' 'text' 'string' 'Add trails before and after regions (s)' } ...
               { 'style' 'edit' 'string' '0.25' } ...
               { 'style' 'text' 'string' 'Use hanning window before computing FFT' } ...
               { 'style' 'checkbox' 'string' '' 'value' 1 } ...
                  };
    geom = { [2 1] [2 1] [2 1] [2 1] [2 1] [2 1] [2 1] };
    result = inputgui('uilist', uilist, 'geometry', geom, 'title', 'Reject continuous portions of data - pop_rejcont()');
    if length( result ) == 0 return; end
    
    options = { 'elecrange'     str2num(result{1}) ...
                'freqlimit'     str2num(result{2}) ...
                'threshold'     str2double(result{3}) ...
                'epochlength'   str2double(result{4}) ...
                'contiguous'    str2double(result{5}) ...
                'addlength'     str2double(result{6}) ...
                'taper'         fastif(result{7}, 'hamming', 'none') };
else
    options = varargin;
end

opt = finputcheck(options, { 'threshold'     { 'real';'cell' }  []    10;
                             'freqlimit'     { 'real';'cell' }  []    [35 128];
                             'mode'          'string' { 'mean';'max' } 'max';
                             'correct'       'string' { 'remove';'blank' } 'remove';
                             'elecrange'     'real'   []    [1:EEG.nbchan];
                             'rejectori'     'real'   []    [];
                             'contiguous'    'real'   []    4;
                             'addlength'     'real'   []    0.25;
                             'precompstruct' 'struct' []    struct([]);
                             'eegplot'       'string' { 'on';'off' } 'off';
                             'onlyreturnselection' 'string' { 'on';'off' } 'off';
                             'verbose'       'string' { 'on';'off' } 'on';
                             'taper'         'string' { 'none' 'hamming' } 'none';
                             'overlap'       'real'   []    0.25;
                             'epochlength'   'real'   []    0.5 }, 'pop_rejcont');
if ischar(opt), error(opt); end
if ~iscell(opt.threshold) && length(opt.threshold) == 2 && ...
    iscell(opt.freqlimit) && length(opt.freqlimit) == 2
    opt.threshold = { opt.threshold(1) opt.threshold(2) };
end
if ~iscell(opt.threshold), opt.threshold = { opt.threshold }; end
if ~iscell(opt.freqlimit), opt.freqlimit = { opt.freqlimit }; end

%EEG.event = [];
grouplen  = opt.contiguous/2*opt.epochlength*EEG.srate+1; % maximum number of points for grouping regions
color     = [ 0 0.9 0];             % color of rejection window

NEWEEG = EEG;
if isempty(opt.precompstruct)
    % compute power spectrum
    % ----------------------
    % average reference 
    % NEWEEG.data(opt.elecrange,:) = NEWEEG.data(opt.elecrange,:)-repmat(mean(NEWEEG.data(opt.elecrange,:),1), [length(opt.elecrange) 1]);

    % only keep boundary events
    % -------------------------
    tmpevent = NEWEEG.event;
    if ~isempty(tmpevent)
        if isnumeric( tmpevent(1).type )
            NEWEEG.event = [];
        else
            boundEvent = strmatch('boundary', { tmpevent.type }, 'exact');
            NEWEEG.event = NEWEEG.event(boundEvent);
        end
    end

    [TMPNEWEEG] = eeg_regepochs(NEWEEG, opt.overlap, [0 opt.epochlength], NaN);
    %[TMPNEWEEG indices] = pop_rejspec(TMPNEWEEG, 1, [1:64], -100, 15, 30, 45, 0, 0);
    %rejepoch = find(indices);
    
    tmpdata = TMPNEWEEG.data;
    if strcmpi(opt.taper, 'hamming'), 
        tmpdata = bsxfun(@times, tmpdata, hamming(size(TMPNEWEEG.data,2))');
    end
    tmp   = fft(tmpdata, [], 2);
    freqs = linspace(0, TMPNEWEEG.srate/2, size(tmp,2)/2);
    freqspectrum = freqs(2:end); % remove DC (match the output of PSD)
    tmp   = tmp(:,2:size(tmp,2)/2,:);
    warning('off', 'MATLAB:log:logOfZero');
    tmpspec  = 10*log10(abs(tmp).^2);  
    warning('on', 'MATLAB:log:logOfZero');
    tmpspec  = tmpspec - repmat( mean(tmpspec,3), [1 1 TMPNEWEEG.trials]);
    specdata = tmpspec;

    % compute mean spectrum
    % ---------------------
    meanspectrum = nan_mean(specdata(opt.elecrange, :, :), 1);
    precompstruct.spec  = meanspectrum;
    precompstruct.freqs = freqspectrum;
    precompstruct.EEG   = TMPNEWEEG;
else
    meanspectrum  = opt.precompstruct.spec;
    freqspectrum  = opt.precompstruct.freqs;
    TMPNEWEEG     = opt.precompstruct.EEG;
    precompstruct = opt.precompstruct;
end

% apply threshold to average of all electrodes
% --------------------------------------------
rejepoch = [];
for iReject = 1:length(opt.threshold)
    threshold = opt.threshold{iReject};
    freqLim   = opt.freqlimit{iReject};
    if length(threshold) == 1, threshold = [ -100 threshold ]; end
    if strcmpi(opt.mode, 'max')
        [I1, tmpRejEpoch, NS, Erej] = eegthresh( meanspectrum, size(meanspectrum,2), 1, threshold(1), threshold(2), [freqspectrum(1) freqspectrum(end)], freqLim(1), freqLim(2));
    else
        tmpRejEpoch = zeros(size(meanspectrum,3),1);
        [~,indb] = min(abs(freqLim(1) - freqspectrum));
        [~,inde] = min(abs(freqLim(2) - freqspectrum));
        for indexe = 1:size(meanspectrum,1)
            sigtmp = meanspectrum(indexe,indb:inde,:);
            if strcmpi(opt.mode, 'mean')
                sigmax = mean(sigtmp,2);
                sigmin = sigmax;
            else
                sigmax = max(sigtmp,2);
                sigmin = min(sigtmp,2);
            end
            sigmin = squeeze(sigmin); % 1 dim at this point
            sigmax = squeeze(sigmax); % 1 dim at this point
            tmpRejEpoch = tmpRejEpoch | ( sigmin < threshold(1) ) | ( sigmax > threshold(2) );
        end
        tmpRejEpoch = find(tmpRejEpoch);
    end
    rejepoch = union_bc(rejepoch, tmpRejEpoch);
    if strcmpi(opt.verbose, 'on')
        fprintf('%d regions selected for rejection, threshold %3.2f-%3.2f dB, frequency limits %3.1f-%3.1f\n', length(tmpRejEpoch), threshold(1), threshold(2), freqLim(1), freqLim(2));
    end
end

% build the winrej array for eegplot
% ----------------------------------
winrej   = [];
if ~isempty(find(cellfun(@isempty, { TMPNEWEEG.event.epoch }) == 1))
    error('Some events are not associated with any epoch');
end
tmpevent = TMPNEWEEG.event;
allepoch = [ tmpevent.epoch ];
if ~isempty(rejepoch)
    for index = 1:length(rejepoch)
        eventepoch = find( rejepoch(index) == allepoch );
        if strcmpi(TMPNEWEEG.event(eventepoch(1)).type, 'X')
            urevent = TMPNEWEEG.event(eventepoch(1)).urevent;
            lat     = TMPNEWEEG.urevent(urevent).latency;
            winrej  = [ winrej; lat lat+opt.epochlength*TMPNEWEEG.srate-1 color ]; %Erej(:,index)'];
        else
            error('Wrong type for epoch');
        end
    end
    winrej(:,6:6+length(opt.elecrange)-1) = 0;
end

% remove isolated regions and merge others
% ----------------------------------------
merged = 0;
isolated = 0;
if opt.contiguous < 2
    for index = size(winrej,1):-1:1
        if index < size(winrej,1) && size(winrej,1) > 1 && winrej(index+1,1) - winrej(index,2) <= grouplen
            winrej(index,2) = winrej(index+1,2);
            winrej(index+1,:) = [];
            merged = merged + 1;
        end
    end
else
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
        end
    end
end
if strcmpi(opt.verbose, 'on')
    fprintf('%d regions merged\n', merged);
    fprintf('%d regions removed\n', isolated);
end

% add time before and after each region
% -------------------------------------
for index = 1:size(winrej,1)
    winrej(index,1) = max(1,         winrej(index,1)-opt.addlength*EEG.srate);
    winrej(index,2) = min(EEG.pnts,  winrej(index,2)+opt.addlength*EEG.srate);
end

% plot result
% -----------
if ~isempty(winrej) 
    selectedregions = winrej(:,1:2);
    if strcmpi(opt.onlyreturnselection, 'off')
        % merge with initial regions
        if ~isempty(opt.rejectori)
            winrej(:,3) = 1; % color
            for iRow = 1:size(opt.rejectori,1)
                winrej(end+1,1:2) = opt.rejectori(iRow,:);
                winrej(end  ,4)   = 1; % color
                winrej(end  ,5)   = 1; % color
            end
        end
        
        command = '[EEG LASTCOM] = pop_select(EEG, ''nopoint'', TMPREJ(:,1:2)); eegh(LASTCOM); [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, ''study'', ~isempty(STUDY)+0); eegh(LASTCOM); eeglab redraw';
        if nargin < 2 || strcmpi(opt.eegplot, 'on')
            eegplot(NEWEEG.data(opt.elecrange,:), 'srate', NEWEEG.srate, 'winrej', winrej, 'command', command, 'events', EEG.event, 'winlength', 50);
            disp('Green is overlap');
            disp('Light blue is ORIGINAL rejection');
            disp('Yellow is AUTOMATIC rejection');
        else
            if strcmpi(opt.correct, 'remove')
                EEG = pop_select(EEG, 'nopoint', round(selectedregions));
            else
                tmpRegions = round(selectedregions);
                for iRegion = 1:size(tmpRegions,1)
                    EEG.data(opt.elecrange, tmpRegions(iRegion,1):tmpRegions(iRegion,2)) = 0;
                end
            end
        end
    else
        EEG = [];
    end
else
    selectedregions = [];
    if strcmpi(opt.verbose, 'on')
        disp('No region removed');
    end
end
if nargout > 3
    com = sprintf('%% the command below does the automated rejection (the pop_select command after that includes manual editing)\n%% EEG = pop_rejcont(EEG, %s);', vararg2str(options));
end
