% std_pac() - Compute or read PAC data (Phase Amplitude Coupling).
%
% Usage:  
%   >> [X times logfreqs ] = std_pac(EEG, 'key', 'val', ...);
% Inputs:
%   EEG          - an EEG dataset structure. 
%
% Optional inputs:
%   'components1'- [numeric vector] components in the EEG structure used 
%                  for spectral amplitude in PAC {default|[]: all }
%   'components2'- [numeric vector] components in the EEG structure used
%                  for phase in PAC {default|[]: all }
%   'channels1'  - [numeric vector or cell array of channel labels] channels 
%                  in the EEG structure for spectral amplitude in PAC 
%                  {default|[]: no channels}
%   'channels2'  - [numeric vector or cell array of channel labels] channels 
%                  in the EEG structure for phase in PAC 
%                  {default|[]: no channels}
%   'freqs'      - [minHz maxHz] the PAC frequency range to compute power.
%                  {default: 12 to EEG sampling rate divided by 2}
%   'cycles'     - [wavecycles (factor)]. If 0 -> DFT (constant window length 
%                  across frequencies).
%                  If >0 -> the number of cycles in each analysis wavelet. 
%                  If [wavecycles factor], wavelet cycles increase with 
%                  frequency, beginning at wavecyles. (0 < factor < 1) 
%                  factor = 0 -> fixed epoch length (DFT, as in FFT). 
%                  factor = 1 -> no increase (standard wavelets)
%                  {default: [0]}
%   'freqphase'  - [valHz] single number for computing the phase at a given
%                  frequency.
%   'cyclephase' - [valcycle] single cycle number.
%   'timewindow' - [minms maxms] time window (in ms) to plot.
%                  {default: all output latencies}
%   'padratio'   - (power of 2). Multiply the number of output frequencies 
%                  by dividing their frequency spacing through 0-padding.
%                  Output frequency spacing is (low_freq/padratio).
%   'recompute'  - ['on'|'off'] 'on' forces recomputation of PAC. 
%                  {default: 'off'}
%
% Other optional inputs:
%   This function will take any of the newtimef() optional inputs (for instance
%   to compute log-space frequencies)...
%
% Outputs:
%   X         - the PAC of the requested ICA components/channels 
%               in the selected frequency and time range. 
%   times     - vector of time points for which the PAC were computed. 
%   freqs     - vector of frequencies (in Hz) at which the 
%               PAC was evaluated. 
%
% Files written or modified:     
%              [dataset_filename].icapac   <-- saved component PAC
% OR for channels
%              [dataset_filename].datpac   <-- saved channel PAC
%
% See also: timef(), std_itc(), std_erp(), std_spec(), std_topo(), std_preclust()
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, July, 2009-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, arno@sccn.ucsd.edu
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

function [X, times, freqs, parameters] = std_pac(EEG, varargin)

if nargin < 1
    help std_pac;
    return;
end

options = {};
[g timefargs] = finputcheck(varargin, { ...
                        'components1'    'integer'     []      [];
                        'channels1'      { 'cell','integer' }  { [],[] }     {};
                        'components2'    'integer'     []      [];
                        'channels2'      { 'cell','integer' }  { [],[] }     {};
                        'outputfile'    'string'      []      '';
                        'powbase'       'real'        []      [];
                        'plot'          'string'      { 'on','off' }      'off';
                        'recompute'     'string'      { 'on','off' }      'off';
                        'getparams'     'string'      { 'on','off' }      'off';
                        'timerange'     'real'        []      []; 
                        'freqrange'     'real'        []      [];
                        'padratio'      'real'        []      1;
                        'freqs'         'real'        []      [12 EEG.srate/2];
                        'cycles'        'real'        []      [8];
                        'freqphase'     'real'        []      [5];
                        'cyclephase'    'real'        []      [3];
                        'interp'        'struct'      { }     struct([]);
                        'rmcomps'       'integer'     []      [];
                        'freqscale'     'string'      []      'log' }, 'std_pac', 'ignore');
if ischar(g), error(g); end

% checking input parameters
% -------------------------
if isempty(g.components1) && isempty(g.channels1)
    if isempty(EEG(1).icaweights)
        error('EEG.icaweights not found');
    end
    g.components1 = 1:size(EEG(1).icaweights,1);
    g.components2 = 1:size(EEG(1).icaweights,1);
    disp('Computing PAC with default values for all components of the dataset');
end

% select ICA components or data channels
% --------------------------------------
if ~isempty(g.outputfile)
    filenamepac   = fullfile('', [ g.outputfile '.datpac' ]);
    g.indices1 = std_chaninds(EEG, g.channels1);
    g.indices2 = std_chaninds(EEG, g.channels2);
    prefix = 'chan';
elseif ~isempty(g.components1)
    g.indices1 = g.components1;
    g.indices2 = g.components2;
    prefix = 'comp';
    filenamepac   = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icapac' ]);
    if ~isempty(g.channels1)
        error('Cannot compute PAC for components and channels at the same time');
    end
elseif ~isempty(g.channels1)
    g.indices1 = std_chaninds(EEG, g.channels1);
    g.indices2 = std_chaninds(EEG, g.channels2);
    prefix = 'chan';
    filenamepac   = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'datpac' ]);
end

% Compute PAC parameters
% -----------------------
parameters = { 'wavelet', g.cycles, 'padratio', g.padratio, ...
               'freqs2', g.freqphase, 'wavelet2', g.cyclephase, 'freqscale', g.freqscale, timefargs{:} };
if length(g.freqs)>0, parameters = { parameters{:} 'freqs' g.freqs }; end

% Check if PAC information found in datasets and if fits requested parameters 
% ----------------------------------------------------------------------------
if exist( filenamepac ) && strcmpi(g.recompute, 'off')
    fprintf('Use existing file for PAC: %s\n', filenamepac);
    if ~isempty(g.components1)
         [X, times, freqs, parameters] = std_readpac(EEG, 1,  g.indices1,  g.indices2, g.timerange, g.freqrange);
    else [X, times, freqs, parameters] = std_readpac(EEG, 1, -g.indices1, -g.indices2, g.timerange, g.freqrange);
    end
    return;
end

% return parameters
% -----------------
if strcmpi(g.getparams, 'on')
    X = []; times = []; freqs = [];
    return;
end

options = {};
if ~isempty(g.components1)
    tmpdata = eeg_getdatact(EEG, 'component', [1:size(EEG(1).icaweights,1)]);
else
    EEG.data = eeg_getdatact(EEG, 'channel', [1:EEG.nbchan], 'rmcomps', g.rmcomps);
    if ~isempty(g.rmcomps), options = { options{:} 'rmcomps' g.rmcomps }; end
    if ~isempty(g.interp), 
        EEG = eeg_interp(EEG, g.interp, 'spherical'); 
        options = { options{:} 'interp' g.interp };
    end
    tmpdata = EEG.data;
end;        

% Compute PAC
% -----------
all_pac = [];
for k = 1:length(g.indices1)  % for each (specified) component
    for l = 1:length(g.indices2)  % for each (specified) component
        tmpparams = parameters;

        % Run pac() to get PAC
        % --------------------
        timefdata1  = tmpdata(g.indices1(k),:,:);
        timefdata2  = tmpdata(g.indices2(l),:,:);
        if strcmpi(g.plot, 'on'), figure; end
        %[logersp,logitc,logbase,times,logfreqs,logeboot,logiboot,alltfX] ...
        [pacvals, times, freqs1, freqs2] = pac( timefdata1, timefdata2, EEG(1).srate, 'tlimits', [EEG.xmin EEG.xmax]*1000, tmpparams{1:end});
        all_pac = setfield( all_pac, [ prefix int2str(g.indices1(k)) '_' int2str(g.indices2(l)) '_pac' ], squeeze(single(pacvals )));
    end
end

% Save PAC into file
% ------------------
all_pac.freqs      = freqs1;
all_pac.times      = times;
all_pac.datatype   = 'PAC';
all_pac.parameters = tmpparams;

if ~isempty(g.channels1)
    if ~isempty(EEG(1).chanlocs)
        tmpchanlocs = EEG(1).chanlocs;
        all_pac.chanlabels1   = { tmpchanlocs(g.indices1).labels };
        all_pac.chanlabels2   = { tmpchanlocs(g.indices2).labels };
    end
end

std_savedat( filenamepac , all_pac );
if ~isempty(g.components1)
     [X, times, freqs, parameters] = std_readpac(EEG, 1,  g.indices1,  g.indices2, g.timerange, g.freqrange);
else [X, times, freqs, parameters] = std_readpac(EEG, 1, -g.indices1, -g.indices2, g.timerange, g.freqrange);
end

% --------------------------------------------------------
% -------------------- READ PAC DATA ---------------------
% --------------------------------------------------------
function [pacvals, freqs, timevals, params] = std_readpac(ALLEEG, abset, comp1, comp2, timewindow, freqrange);

if nargin < 5
    timewindow = [];
end
if nargin < 6
    freqrange = [];
end

% multiple entry
% --------------
if length(comp1) > 1 || length(comp2) > 1
    for index1 = 1:length(comp1)
        for index2 = 1:length(comp2)
            [tmppac, freqs, timevals, params] = std_readpac(ALLEEG, abset, comp1(index1), comp2(index2), timewindow, freqrange);
            pacvals(index1,index2,:,:,:) = tmppac;
        end
    end
    return;
end

for k = 1: length(abset)    
    
    if comp1 < 0
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'datpac']);
        comp1   = -comp1;
        comp2   = -comp2;
        prefix = 'chan';
    else    
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'icapac']);
        prefix = 'comp';
    end
    try
        tmppac = load( '-mat', filename, 'parameters', 'times', 'freqs');
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end
    
    tmppac.parameters = removedup(tmppac.parameters);
    params    = struct(tmppac.parameters{:});
    params.times = tmppac.times;
    params.freqs = tmppac.freqs;
    if isempty(comp1)
        pacvals    = [];
        freqs  = [];
        timevals  = [];
        return;
    end
    tmppac   = load( '-mat', filename, 'parameters', 'times', 'freqs', ...
                     [ prefix int2str(comp1) '_' int2str(comp2) '_pac']);
    
    pacall{k} = double(getfield(tmppac, [ prefix int2str(comp1) '_' int2str(comp2) '_pac']));
    tlen      = length(tmppac.times);
    flen      = length(tmppac.freqs);

end

% select plotting or clustering time/freq range
% ---------------------------------------------
if ~isempty(timewindow)
    if timewindow(1) > tmppac.times(1) || timewindow(end) < tmppac.times(end)
        maxind = max(find(tmppac.times <= timewindow(end)));
        minind = min(find(tmppac.times >= timewindow(1)));
    else
        minind = 1;
        maxind = tlen;
    end
else
    minind = 1;
    maxind = tlen;
end
if ~isempty(freqrange)
    if freqrange(1) > exp(1)^tmppac.freqs(1) || freqrange(end) < exp(1)^tmppac.freqs(end)
        fmaxind = max(find(tmppac.freqs <= freqrange(end)));
        fminind = min(find(tmppac.freqs >= freqrange(1)));
    else
        fminind = 1;
        fmaxind = flen;
    end
else
    fminind = 1;
    fmaxind = flen;
end

% return parameters
% ----------------
for cond  = 1:length(abset)
    try
        pac = pacall{cond}(fminind:fmaxind,minind:maxind);
    catch
        pac = pacall{cond}; % for 'method', 'latphase'
    end
    pacvals(:,:,cond) = pac;
end
freqs    = tmppac.freqs(fminind:fmaxind);
timevals = tmppac.times(minind:maxind);

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
    [tmp indices] = unique_bc(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        %fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end
    cella = cella(sort(union(indices*2-1, indices*2)));

