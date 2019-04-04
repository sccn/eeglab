% std_specgram() - Returns the ICA component or channel spectrogram for a dataset.
%                  Saves the spectra in a file.
% Usage:    
%           >> [spec freqs] = std_specgram(EEG, 'key', 'val', ...);
%
% Inputs:
%   EEG - a loaded epoched EEG dataset structure. 
%
% Optional inputs:
%   'components' - [numeric vector] components of the EEG structure for which 
%                  activation spectogram will be computed. Note that because 
%                  computation of component spectra is relatively fast, all 
%                  components spectra are computed and saved. Only selected 
%                  component are returned by the function to Matlab
%                  {default|[] -> all}
%   'channels'   - [cell array] channels of the EEG structure for which 
%                  activation spectogram will be computed. Note that because 
%                  computation of spectrum is relatively fast, all channels 
%                  spectrum are computed and saved. Only selected channels 
%                  are returned by the function to Matlab
%                  {default|[] -> none}
%   'recompute'  - ['on'|'off'] force recomputing ERP file even if it is 
%                  already on disk.
%
% Other optional spectral parameters:
%   All optional parameters to the newtimef function may be provided to this function
%   as well.
%
% Outputs:
%   spec      - the mean spectra (in dB) of the requested ICA components in the selected 
%               frequency range (with the mean of each spectrum removed). 
%   freqs     - a vector of frequencies at which the spectra have been computed. 
%
% Files output or overwritten for ICA: 
%               [dataset_filename].icaspecgram,
% Files output or overwritten for data: 
%               [dataset_filename].datspecgram, 
% 
% See also  spectopo(), std_erp(), std_ersp(), std_map(), std_preclust()
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, January, 2005

% Defunct:      0 -> if frequency range is different from saved spectra, ask via a 
%                    pop-up window whether to keep existing spectra or to overwrite them. 

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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

% eeg_specgram() - Compute spectrogramme taking into account boundaries in
%                  the data.
% Usage:
%   >> EEGOUT = eeg_specgram( EEG, typeplot, num, 'key', 'val');
%
% Inputs:
%   EEG      - EEG dataset structure
%   typeplot - type of processinopt. 1 process the raw
%              data and 0 the ICA components
%   num      - component or channel number
%
% Optional inputs:
%   'winsize' - [integer] window size in points
%   'overlap' - [integer] window overlap in points (default: 0)
%   'movav'   - [real] moving average
%
% Author: Arnaud Delorme, CERCO, CNRS, 2008-

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [erspinterp t f ] = eeg_specgram(EEG, varargin);

if nargin < 1
    help std_specgram;
    return;
end

[opt moreopts] = finputcheck(varargin, { 'components' 'integer' []             [];
                                         'channels'   { 'cell','integer' }  { [] [] }     {}
                                         'recompute'  'string'  { 'on','off' } 'off';
                                         'winsize'    'integer' []             3; % 3 seconds
                                         'rmcomps'    'integer' []             [];
                                         'interp'     'struct'  { }            struct([]);
                                         'overlap'    'integer' []             0;  
                                         'plot'       'string'  { 'off','on' } 'off';
                                         'freqrange'  'real'    []             [];
                                         'timerange'  'real'    []             [];
                                         'filter'     'real'    []             []}, ...    % 11 points
                                         'eeg_specgram', 'ignore');
if ischar(opt), error(opt); end
if isfield(EEG,'icaweights')
   numc = size(EEG.icaweights,1);
else
   error('EEG.icaweights not found');
end
if isempty(opt.components)
    opt.components = 1:numc;
end
%opt.winsize = 2^ceil(log2(opt.winsize*EEG.srate));
opt.winsize = opt.winsize*EEG.srate;

% filename 
% --------
if ~isempty(opt.channels)
    filename = fullfile( EEG.filepath,[ EEG.filename(1:end-3) 'datspecgram']);
    prefix = 'chan';
    opt.indices = opt.channels;
    if iscell(opt.channels)
        tmpchanlocs = EEG(1).chanlocs;
        for index = 1:length(opt.channels)
            chanind = strmatch( lower(opt.channels{index}), lower({ tmpchanlocs.labels }), 'exact');
            if isempty(chanind), error('Channel group not found'); end
            chaninds(index) = chanind;
        end
        opt.indices  = chaninds;
        opt.channels = chaninds;
    end
else    
    filename = fullfile( EEG.filepath,[ EEG.filename(1:end-3) 'icaspecgram']);
    prefix = 'comp';
    opt.indices = opt.components;
end

% SPEC information found in datasets
% ----------------------------------
if exist(filename) && strcmpi(opt.recompute, 'off')

    if strcmpi(prefix, 'comp')
        [erspinterp, t, f] = std_readspecgram(EEG, 1, opt.components, opt.freqrange);
    else
        [erspinterp, t, f] = std_readspecgram(EEG, 1, -opt.channels, opt.freqrange);
    end
    return;
    
end

% No SPEC information found
% ------------------------
options = {};
if strcmpi(prefix, 'comp')
    X = eeg_getdatact(EEG, 'component', [1:size(EEG.icaweights,1)]);
else
    EEG.data = eeg_getdatact(EEG, 'channel', [1:EEG.nbchan], 'rmcomps', opt.rmcomps);
    if ~isempty(opt.rmcomps), options = { options{:} 'rmcomps' opt.rmcomps }; end
    if ~isempty(opt.interp), 
        EEG = eeg_interp(EEG, opt.interp, 'spherical'); 
    end
    X = EEG.data;
end

% get the array of original point latency
% ---------------------------------------
urpnts      = eeg_urpnts(EEG);
urarray     = eeg_makeurarray(EEG, urpnts); % contain the indices of the urpoint in the EEG data
                                            % urarray(1000) = 1000, urarray(2300) = 1600 if part removed in the data
urwincenter = opt.winsize/2+1:opt.winsize-opt.overlap:urpnts-opt.winsize/2;
wintag      = ones(1, length(urwincenter));
if EEG.trials == 1
    for i = 1:length(urwincenter)
        win = urwincenter(i)+[-opt.winsize/2+1:opt.winsize/2];
        if ~all(urarray(win))
            wintag(i) = 0;
            %fprintf('Missing data window: %3.1f-%3.1f s\n', (win(1)-1)/EEG.srate, (win(end)-1)/EEG.srate); 
        end
    end
else
    error('eeg_specgram can only be run on continuous data');
end

% compute spectrum 2 solutions
% 1- use newtimef, have to set the exact times and window
% 2- redo the FFT myself
% ----------------------
wincenter = urwincenter(find(wintag));       % remove bad windows
wincenter = urarray(wincenter);              % latency in current dataset
wincenter = 1000*(wincenter-1)/EEG.srate;    % convert to ms
freqs     = linspace(0.1, 50, 100);
options   = { 0 'winsize', opt.winsize, 'baseline', [0 Inf], 'timesout', wincenter, ...
              'plotersp', 'off', 'plotitc', 'off', 'freqs', freqs };
%freqs     = exp(linspace(log(EEG.srate/opt.winsize*4), log(50), 100));
%cycles    = linspace(3,8,100);
%options   = { [3 0.8] 'winsize', opt.winsize, 'baseline', [0 Inf], 'timesout', wincenter, ...
%              'freqs' freqs 'cycles' cycles 'plotersp', 'off', 'plotitc', 'off' };
for ic = 1:length(opt.indices)
    [ersp(:,:,ic) itc powebase t f] = newtimef(X(opt.indices(ic), :), EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, options{:}, moreopts{:});
end

% interpolate and smooth in time
% ------------------------------
disp('Now interpolating...');
wininterp  = find(wintag == 0);
erspinterp = zeros(size(ersp,1), length(urwincenter), size(ersp,3));
erspinterp(:,find(wintag),:) = ersp;
for s = 1:size(ersp,3)
    for i=1:length(wininterp)
        first1right = find(wintag(wininterp(i):end));
        first1left  = find(wintag(wininterp(i):-1:1));
        if isempty(first1right)
            erspinterp(:,wininterp(i),s) = erspinterp(:,wininterp(i)+1-first1left(1),s);
        elseif isempty(first1left)
            erspinterp(:,wininterp(i),s) = erspinterp(:,wininterp(i)-1+first1right(1),s);
        else
            erspinterp(:,wininterp(i),s) =(erspinterp(:,wininterp(i)-1+first1right(1),s) + erspinterp(:,wininterp(i)+1-first1left(1),s))/2;
        end
    end
end

%erspinterp = vectdata(ersp, urwincenter(find(wintag))/EEG.srate, 'timesout', urwincenter/EEG.srate);

% smooth in time with a simple convolution
% ----------------------------------------
if ~isempty(opt.filter)
    filterlen = opt.filter(1);
    filterstd = opt.filter(2);
    incr = 2*filterstd/(filterlen-1); %gaussian filter
    filter = exp(-(-filterstd:incr:filterstd).^2);
    
    erspinterp = convn(erspinterp, filter/sum(filter), 'same');
    %erspinterp = conv2(erspinterp, filter/sum(filter)); 
    %erspinterp(:, [1:(filterlen-1)/2 end-(filterlen-1)/2+1:end]) = [];        
end

% plot result
% -----------
t = (urwincenter-1)/EEG.srate;
if strcmpi(opt.plot, 'on')
    figure; imagesc(t, log(f), erspinterp);
    ft = str2num(get(gca,'yticklabel'));
    ft = exp(1).^ft;
    ft = unique_bc(round(ft));
    ftick = get(gca,'ytick');
    ftick = exp(1).^ftick;
    ftick = unique_bc(round(ftick));
    ftick = log(ftick);
    set(gca,'ytick',ftick);
    set(gca,'yticklabel', num2str(ft));

    xlabel('Time (h)');
    ylabel('Frequency (Hz)');
    set(gca, 'ydir', 'normal');
end

% Save SPECs in file (all components or channels)
% ----------------------------------
options = { 'winsize' opt.winsize 'overlap' opt.overlap moreopts{:} };
if strcmpi(prefix, 'comp')
    savetofile( filename, t, f, erspinterp, 'comp', opt.indices, options, [], opt.interp);
    [erspinterp, t, f] = std_readspecgram(EEG, 1, opt.components, opt.timerange, opt.freqrange);
else
    tmpchanlocs = EEG(1).chanlocs;
    savetofile( filename, t, f, erspinterp, 'chan', opt.indices, options, { tmpchanlocs.labels }, opt.interp);
    [erspinterp, t, f] = std_readspecgram(EEG, 1, -opt.channels, opt.timerange, opt.freqrange);
end
return;

% recompute the original data length in points
% --------------------------------------------
function urlat = eeg_makeurarray(EEG, urpnts);
    
    if isempty(EEG.event) || ~isfield(EEG.event, 'duration')
        urlat = 1:EEG.pnts;
        return;
    end
    
    % get boundary events latency and duration
    % ----------------------------------------
    tmpevent  = EEG.event;
    bounds    = strmatch('boundary', { tmpevent.type });
    alldurs   = [ tmpevent(bounds).duration ];
    alllats   = [ tmpevent(bounds).latency  ];
    if length(alldurs) >= 1
        if alldurs(1) <= 1
            alllats(1) = [];
            alldurs(1) = [];
        end
    end
    
    if isempty(alllats)
        urlat = 1:EEG.pnts;
        return;
    end
    
    % build the ur boolean array
    % --------------------------
    urlat = ones(1, urpnts);
    for i=1:length(alllats)
        urlat(round(alllats(i)+0.5):round(alllats(i)+0.5+alldurs(i)-1)) = 0;
        alllats(i+1:end) = alllats(i+1:end)+alldurs(i);
    end
    urlat(find(urlat)) = 1:EEG.pnts;

% -------------------------------------
% saving SPEC information to Matlab file
% -------------------------------------
function savetofile(filename, t, f, X, prefix, comps, params, labels, interp);
    
    disp([ 'Saving SPECTRAL file ''' filename '''' ]);
    allspec = [];
    for k = 1:length(comps)
        allspec = setfield( allspec, [ prefix int2str(comps(k)) ], X(:,:,k));
    end
    if ~isempty(labels)
        allspec.labels = labels;
    end
    allspec.freqs      = f;
    allspec.times      = t;
    allspec.parameters = params;
    allspec.datatype   = 'SPECTROGRAM';
    allerp.interpolation = fastif(isempty(interp), 'no', interp);
    allspec.average_spec = mean(X,1);
    std_savedat(filename, allspec);
    
% recompute the original data length in points
% --------------------------------------------
function pntslat = eeg_urpnts(EEG);
    
    if isempty(EEG.event) || ~isfield(EEG.event, 'duration')
        pntslat = EEG.pnts;
        return;
    end
    tmpevent = EEG.event;
    bounds = strmatch('boundary', { tmpevent.type });
    alldurs = [ tmpevent(bounds).duration ];
    if length(alldurs) > 0
        if alldurs(1) <= 1, alldurs(1) = [];
        end
    end
    pntslat = EEG.pnts + sum(alldurs);

% recompute the original latency
% ------------------------------
function pntslat = eeg_urlatency(EEG, pntslat);

    if isempty(EEG.event), return; end
    if ~ischar(EEG.event(1).type), return; end
    
    tmpevent = EEG.event;
    bounds = strmatch('boundary', { tmpevent.type })
    for i=1:length(bounds)
        if EEG.event(bounds(i)).duration > 1
            pntslat = pntslat + EEG.event(bounds(i)).duration;
        end
    end
    
