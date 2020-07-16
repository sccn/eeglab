% std_ersp() - Compute ERSP and/or ITC transforms for ICA components 
%              or data channels of a dataset. Save results into Matlab 
%              float files. 
%
% Function description:
%              The function computes the mean ERSP or ITC for the selected 
%              dataset ICA components or data channels in the requested 
%              frequency range and time window (the two are dependent). 
%              Frequencies are equally log spaced. Options specify component 
%              numbers, desired frequency range,  time window length, 
%              frequency resolution, significance level, and wavelet
%              cycles. See >> help newtimef and >> timef details 
%
%              Two Matlab files are saved (for ERSP and ITC). These contain 
%              the ERSP|ITC image, plus the transform parameters 
%              used to compute them. Saves the computed dataset mean images 
%              in dataset-name files with extensions '.icaersp' and '.icaitc'
%              for ICA components or '.datersp', '.datitc' for data channels.
% Usage:  
%              >> [X times logfreqs ] = std_ersp(EEG, 'key', 'val', ...);
% Inputs:
%   EEG          - a loaded epoched EEG dataset structure. May be an array
%                  of such structure containing several datasets.
%
% Other inputs:
%   'trialindices' - [cell array] indices of trials for each dataset.
%                  Default is EMPTY (no trials). NEEDS TO BE SET.
%   'components' - [numeric vector] components of the EEG structure for which 
%                  activation spectrum will be computed. Note that because 
%                  computation of ERP is so fast, all components spectrum are
%                  computed and saved. Only selected component 
%                  are returned by the function to Matlab
%                  {default|[] -> all}
%   'channels'   - [cell array] channels of the EEG structure for which 
%                  activation spectrum will be computed. Note that because 
%                  computation of ERP is so fast, all channels spectrum are
%                  computed and saved. Only selected channels 
%                  are returned by the function to Matlab
%                  {default|[] -> none}
%   'recompute'  - ['on'|'off'] force recomputing ERP file even if it is 
%                  already on disk.
%   'recompute'  - ['on'|'off'] force recomputing data file even if it is 
%                  already on disk.
%   'rmcomps'    - [integer array] remove artifactual components (this entry
%                  is ignored when plotting components). This entry contains 
%                  the indices of the components to be removed. Default is none.
%   'interp'     - [struct] channel location structure containing electrode
%                  to interpolate ((this entry is ignored when plotting 
%                  components). Default is no interpolation.
%   'fileout'    - [string] name of the file to save on disk. The default
%                  is the same name (with a different extension) as the 
%                  dataset given as input.
%  'savetrials'  - ['on'|'off'] save single-trials ERSP. Requires a lot of disk
%                  space (dataset space on disk times 10) but allow for refined
%                  single-trial statistics. This option is obsolete as
%                  trials are now always saved.
%  'savefile'    - ['on'|'off'] save file or simply return measures.
%                  Default is to save files ('on').
%  'getparams'   - ['on'|'off'] return optional parameters for the newtimef 
%                  function (and do not compute anything). This argument is
%                  obsolete (default is 'off').
%
% ERSP optional inputs:
%   'type'       - ['ersp'|'itc'|'ersp&itc'] save ERSP, ITC, or both data 
%                  types to disk {default: 'ersp'}
%   'freqs'      - [minHz maxHz] the ERSP/ITC frequency range to compute 
%                  and return. {default: 3 to EEG sampling rate divided by 3}
%   'timelimits' - [minms maxms] time window (in ms) to compute.
%                  {default: whole input epoch}.
%   'cycles'     - [wavecycles (factor)]. If 0 -> DFT (constant window length 
%                  across frequencies).
%                  If >0 -> the number of cycles in each analysis wavelet. 
%                  If [wavecycles factor], wavelet cycles increase with 
%                  frequency, beginning at wavecyles. (0 < factor < 1) 
%                  factor = 0 -> fixed epoch length (DFT, as in FFT). 
%                  factor = 1 -> no increase (standard wavelets)
%                  {default: [0]}
%   'padratio'   - (power of 2). Multiply the number of output frequencies 
%                  by dividing their frequency spacing through 0-padding.
%                  Output frequency spacing is (low_freq/padratio).
%   'alpha'      - If in (0, 1), compute two-tailed permutation-based 
%                  probability thresholds and use these to mask the output 
%                  ERSP/ITC images {default: NaN}
%   'powbase'    - Deprecated. Note that baseline can be readjusted after 
%                  computation as single trial spectral decompositions are stored.
%
% Other optional inputs:
%   This function will take any of the newtimef() optional inputs (for instance
%   to compute log-space frequencies)...
%
% Outputs:
%   X         - the masked log ERSP/ITC of the requested ICA components/channels 
%               in the selected frequency and time range. Note that for
%               optimization reasons, this parameter is now empty or 0. X
%               thus must be read from the datafile saved on disk.
%   times     - vector of time points for which the ERSPs/ITCs were computed. 
%   logfreqs  - vector of (equally log spaced) frequencies (in Hz) at which the 
%               log ERSP/ITC was evaluated. 
%   parameters - parameters given as input to the newtimef function.
%
% Files written or modified:     
%              [dataset_filename].icaersp   <-- saved component ERSPs
%              [dataset_filename].icaitc    <-- saved component ITCs
%              [dataset_filename].icatimef  <-- saved component single
%                                               trial decompositions.
%  OR for channels
%              [dataset_filename].datersp   <-- saved channel ERSPs
%              [dataset_filename].datitc    <-- saved channel ITCs
%              [dataset_filename].dattimef  <-- saved channel single
%                                               trial decompositions.
% Example: 
%            % Create mean ERSP and ITC images on disk for all comps from 
%            % dataset EEG use three-cycle wavelets (at 3 Hz) to more than 
%            % three-cycle wavelets at 50 Hz. See >> help newtimef 
%            % Return the (equally log-freq spaced, probability-masked) ERSP.
%            >> [Xersp, times, logfreqs] = std_ersp(EEG, ...
%                       'type', 'ersp', 'freqs', [3 50], 'cycles', [3 0.5]);
%
% See also: timef(), std_itc(), std_erp(), std_spec(), std_topo(), std_preclust()
%
% Authors: Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, January, 2005-

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

function [X, times, logfreqs, parameters] = std_ersp(EEG, varargin)

if nargin < 1
    help std_ersp;
    return;
end

X = [];
options = {};
if length(varargin) > 1 
    if ~ischar(varargin{1})
        if length(varargin) > 0, options = { options{:} 'components' varargin{1} }; end
        if length(varargin) > 1, options = { options{:} 'freqs'      varargin{2} }; end
        if length(varargin) > 2, options = { options{:} 'timewindow' varargin{3} }; end
        if length(varargin) > 3, options = { options{:} 'cycles'     varargin{4} }; end
        if length(varargin) > 4, options = { options{:} 'padratio'   varargin{5} }; end
        if length(varargin) > 5, options = { options{:} 'alpha'      varargin{6} }; end
        if length(varargin) > 6, options = { options{:} 'type'       varargin{7} }; end
        if length(varargin) > 7, options = { options{:} 'powbase'    varargin{8} }; end
    else
        options = varargin;
    end
end

[g, timefargs] = finputcheck(options, { ...
                        'components'    'integer'               []          [];
                        'channels'      { 'cell','integer' }    { [] [] }   {};
                        'powbase'       'real'                  []          [];
                        'trialindices' { 'integer','cell' }     []          [];
                        'savetrials'    'string'      { 'on','off' }        'off';
                        'plot'          'string'      { 'on','off' }        'off'; % not documented for debugging purpose
                        'recompute'     'string'      { 'on','off' }        'off';
                        'getparams'     'string'      { 'on','off' }        'off';
                        'savefile'      'string'      { 'on','off' }        'on';
                        'parallel'      'string'      { 'on','off' }        'off';
                        'timewindow'    'real'                  []          [];    % ignored, deprecated
                        'fileout'       'string'                []          '';
                        'timelimits'    'real'                  []          [EEG(1).xmin EEG(1).xmax]*1000;
                        'cycles'        'real'                  []          [3 .5];
                        'padratio'      'real'                  []          1;
                        'trialinfo'     'struct'                []          struct([]);
                        'freqs'         'real'                  []          [0 EEG(1).srate/2];
                        'rmcomps'       'cell'                  []          cell(1,length(EEG));
                        'interp'        'struct'                { }         struct([]);
                        'freqscale'     'string'                []         'log';
                        'alpha'         'real'                  []          NaN;
                        'baseline'      'real'                  []          0;
                        'type'          'string'      { 'ersp','itc','both','ersp&itc' }  'both'}, 'std_ersp', 'ignore');
if ischar(g), error(g); end
if isempty(g.trialindices), g.trialindices = cell(length(EEG)); end
if ~iscell(g.trialindices), g.trialindices = { g.trialindices }; end
if ~isempty(g.powbase), disp('''powbase'' parameter is no longer supported at computation time'); end

% checking input parameters
% -------------------------
if isempty(g.components) && isempty(g.channels)
    if isempty(EEG(1).icaweights)
        error('EEG.icaweights not found');
    end
    g.components = 1:size(EEG(1).icaweights,1);
    disp('Computing ERSP with default values for all components of the dataset');
end

% select ICA components or data channels
% --------------------------------------
if isempty(g.fileout), g.fileout = fullfile(EEG(1).filepath, EEG(1).filename(1:end-4)); end
if ~isempty(g.components)
    g.indices = g.components;
    prefix = 'comp';
    filenameersp   = [ g.fileout '.icaersp'  ];
    filenameitc    = [ g.fileout '.icaitc'   ];
    filenametrials = [ g.fileout '.icatimef' ];    
    if ~isempty(g.channels)
        error('Cannot compute ERSP/ITC for components and channels at the same time');
    end
elseif ~isempty(g.channels)
    if iscell(g.channels)
        if ~isempty(g.interp)
            g.indices = eeg_chaninds(g.interp, g.channels, 0);
        else
            g.indices = eeg_chaninds(EEG(1), g.channels, 0);
            for ind = 2:length(EEG)
                if ~isequal(eeg_chaninds(EEG(ind), g.channels, 0), g.indices)
                    error([ 'Channel information must be consistant when ' 10 'several datasets are merged for a specific design' ]);
                end
            end
        end
    else
        g.indices = g.channels;
    end
    prefix = 'chan';
    filenameersp   = [ g.fileout '.datersp'  ];
    filenameitc    = [ g.fileout '.datitc'   ];
    filenametrials = [ g.fileout '.dattimef' ];    
end

% Check if ERSP/ITC information found in datasets and if fits requested parameters 
% ----------------------------------------------------------------------------
if exist( filenameersp ) && strcmpi(g.recompute, 'off')
    fprintf('Use existing file for ERSP: %s; check the ''recompute checkbox'' to force recomputing.\n', filenameersp);
    return;
end

% Compute ERSP parameters
% -----------------------
parameters = { 'cycles', g.cycles, 'padratio', g.padratio, ...
               'alpha', g.alpha, 'freqscale', g.freqscale, timefargs{:} };
defaultlowfreq = 3;
[time_range] = compute_ersp_times(g.cycles,  EEG(1).srate, ...
                                 [EEG(1).xmin EEG(1).xmax]*1000 , defaultlowfreq, g.padratio); 
if time_range(1) < time_range(2) && g.freqs(1) == 0
     g.freqs(1) = defaultlowfreq; % for backward compatibility
end
parameters = { parameters{:} 'freqs' g.freqs };
if strcmpi(g.plot, 'off')
    parameters = { parameters{:} 'plotersp', 'off', 'plotitc', 'off', 'plotphase', 'off' };
end
parameters{end+1} = 'baseline';
parameters{end+1} = g.baseline;

% return parameters
% -----------------
if strcmpi(g.getparams, 'on')
    X = []; times = []; logfreqs = [];
    if strcmpi(g.savetrials, 'on')
        parameters = { parameters{:} 'savetrials', g.savetrials };
    end
    return;
end

options = {};
if ~isempty(g.rmcomps), options = { options{:} 'rmcomps' g.rmcomps }; end
if ~isempty(g.interp),  options = { options{:} 'interp' g.interp }; end
if isempty(g.channels)
     X = eeg_getdatact(EEG, 'component', g.indices, 'trialindices', g.trialindices );
else X = eeg_getdatact(EEG, 'channel'  , g.indices, 'trialindices', g.trialindices, 'rmcomps', g.rmcomps, 'interp', g.interp);
end
if size(X, 3) == 1
    error('The data is continuous for one of the dataset. ERSP can only be computed when data trials are present');
end

% frame range
% -----------
pointrange1 = round(max((g.timelimits(1)/1000-EEG(1).xmin)*EEG(1).srate, 1));
pointrange2 = round(min(((g.timelimits(2)+1000/EEG(1).srate)/1000-EEG(1).xmin)*EEG(1).srate, EEG(1).pnts));
pointrange = [pointrange1:pointrange2];

% Compute ERSP && ITC
% ------------------
allTrialsTmp   = cell(1,length(g.indices));
allTrialsTime  = cell(1,length(g.indices));
allTrialsFreqs = cell(1,length(g.indices));
eeglab_options;
usesingle = option_single;

disp('Computing time/frequency decomposition...');

% CHANGE THE LINE BELOW TO PARFOR TO USE THE PARALLEL TOOLBOX
% numWorkers = 0;
% try
%     if strcmpi(g.parallel, 'on') && ~exist('parpool') && exist('gcp')
%         poolobj = gcp('nocreate');
%         if ~isempty(poolobj)
%             numWorkers = poolobj.NumWorkers;
%         end
%     elseif strcmpi(g.parallel, 'on')
%         try
%             myCluster = parcluster('local');
%             numWorkers = myCluster.NumWorkers;
%         catch, disp('Cound not start parallel job; ERSP will still be computed'); end
%     end
% catch
%     g.parallel = 'off';
% end
%parfor (k = 1:length(g.indices),numWorkers)  % for each (specified) component/channel
for k = 1:length(g.indices)
    tmpparams = parameters;
    if length(g.indices) > 1
        tmpparams{end+1} = 'verbose';
        tmpparams{end+1} = 'off';
    end
    
    % Run timef() to get ERSP
    % ------------------------
    timefdata  = reshape(X(k,pointrange,:), 1, length(pointrange)*size(X,3));
    mytimes = [];
    mylogfreqs = [];
    alltfX   = [];
    if ~isempty(timefdata)
        [logersp,logitc,logbase,mytimes,mylogfreqs,logeboot,logiboot,alltfX] ...
              = newtimef( timefdata, length(pointrange), g.timelimits, EEG(1).srate, tmpparams{2:end});
        %figure; newtimef( TMP.data(32,:), EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles, 'freqs', freqs);
        %figure; newtimef( timefdata, length(pointrange), g.timelimits, EEG.srate, cycles, 'freqs', freqs);
    end
    %if strcmpi(g.plot, 'on'), return; end
    if usesingle
        alltfX = single(alltfX);
    end

    allTrialsTmp{k}   = single( alltfX );
    allTrialsTime{k}  = mytimes;
    allTrialsFreqs{k} = mylogfreqs;
end
all_trials = [];
for k = 1:length(g.indices)  % for each (specified) component/channel
    all_trials = setfield( all_trials, [ prefix int2str(g.indices(k)) ], allTrialsTmp{k});
end
X = allTrialsTmp{1};

% Save ERSP into file
% -------------------
logfreqs             = allTrialsFreqs{1};
times                = allTrialsTime{1};
all_trials.freqs     = allTrialsFreqs{1};
all_trials.times     = allTrialsTime{1};
all_trials.parameters = { options{:} parameters{:} };
all_trials.datatype   = 'TIMEF';
all_trials.datafiles  = computeFullFileName( { EEG.filepath }, { EEG.filename });
all_trials.datatrials = g.trialindices;

all_trials.parameters = parameters;
if ~isempty(g.channels)
    if ~isempty(g.interp)
        all_trials.labels = { g.interp(g.indices).labels };
    elseif ~isempty(EEG(1).chanlocs)
        tmpchanlocs = EEG(1).chanlocs;
        all_trials.labels = { tmpchanlocs(g.indices).labels };
    end
end
all_trials.trialinfo = g.trialinfo;

if strcmpi(g.savefile, 'on')
    std_savedat( filenametrials , all_trials );
end

% compute full file names
% -----------------------
function res = computeFullFileName(filePaths, fileNames);
for index = 1:length(fileNames)
    res{index} = fullfile(filePaths{index}, fileNames{index});
end

