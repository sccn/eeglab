% std_ersp() - Compute ERSP and/or ITC transforms for ICA components 
%              or data channels of a dataset. Save results into Matlab 
%              float files. When these output files already exist, loads 
%              the ERSP/ITC information from them unless the requested 
%              flag specifies differently. If so, a query window 
%              pops up.
%
% Function description:
%              The function returns the masked (as per the requested alpha) 
%              mean ERSP or ITC for the selected dataset ICA components or 
%              data channels in the requested frequency range and time window 
%              (the two are dependent). Frequencies are equally log spaced.
%              Options specify component numbers, desired frequency range, 
%              time window length, frequency resolution, significance level, 
%              and wavelet cycles. See >> help newtimef and >> timef details 
%
%              Two Matlab files are saved (for ERSP and ITC). These contain 
%              the ERSP|ITC image, plus the transform parameters 
%              used to compute them. Saves the computed dataset mean images 
%              in dataset-name files with extensions '.icaersp' and '.icaitc'
%              for ICA components or '.datersp', '.datitc' for data channels.
%              If the ERSPs/ITCs were previously saved into these files, 
%              and the same set of ERSP/ITC parameters are used, the values 
%              are not recomputed, but the information is read from these
%              files. Vectors of frequencies and latencies for the ERSP/ITC 
%              images are returned separately. Returned 'EEG.etc' fields
%              are modified with pointers to the output float files and some 
%              information about them. 
% Usage:  
%              >> [X times logfreqs ] = std_ersp(EEG, 'key', 'val', ...);
% Inputs:
%   EEG          - an EEG dataset structure. 
%
% Optional inputs:
%   'components' - [numeric vector] components in the EEG structure for which 
%                  ERSP and ITC data will be computed {default|[]: all 
%                  components if no 'channels' are specified (see below)}
%   'channels'   - [numeric vector or cell array of channel labels] channels 
%                  in the EEG structure for which ERSP and ITC will be computed 
%                  {default|[]: no channels}
%   'freqs'      - [minHz maxHz] the ERSP/ITC frequency range to compute 
%                  and return. {default: 3 to EEG sampling rate divided by 3}
%   'timelimits' - [minms maxms] time window (in ms) to compute.
%                  {default: whole input epoch}.
%   'timewindow' - [minms maxms] time window (in ms) to plot.
%                  {default: all output latencies}
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
%   'type'       - ['ersp'|'itc'] though both ERSP and ITC images are computed 
%                  and saved to disk, only this transform is returned to the 
%                  command line (see first output, X, below) {default: 'ersp'}
%   'savetrials' - ['on'|'off'] Save single-trial time-freq. decompositions in
%                  a file with extension '.dattimef' (channels) or '.icatimef' 
%                  (components). {default: 'off'}
%   'powbase'    - [ncomps,nfreqs] optional input matrix giving baseline power 
%                  spectra (not dB power, see >> help timef). 
%                  For use in repeated calls to timef() using the same baseine
%                  {default|[] -> none; data windows centered before 0 latency}
%   'recompute'  - ['on'|'off'] 'on' forces recomputation of both ERSP and ITC. 
%                  {default: 'off'}
%
% Other optional inputs:
%   This function will take any of the newtimef() optional inputs (for instance
%   to compute log-space frequencies)...
%
% Outputs:
%   X         - the masked log ERSP/ITC of the requested ICA components/channels 
%               in the selected frequency and time range. 
%   times     - vector of time points for which the ERSPs/ITCs were computed. 
%   logfreqs  - vector of (equally log spaced) frequencies (in Hz) at which the 
%               log ERSP/ITC was evaluated. 
%
% Files written or modified:     
%              [dataset_filename].icaersp   <-- saved component ERSPs
%              [dataset_filename].icaitc    <-- saved component ITCs
%  OR for channels
%              [dataset_filename].datersp   <-- saved channel ERSPs
%              [dataset_filename].datitc    <-- saved channel ITCs
% Example: 
%            % Create mean ERSP and ITC images on disk for all comps from 
%            % dataset EEG use three-cycle wavelets (at 3 Hz) to more than 
%            % three-cycle wavelets at 50 Hz. Use probability masking at 
%            % p < 0.01 with padratio 4. See >> timef details. 
%            % Return the (equally log-freq spaced, probability-masked) ERSP.
%            >> [Xersp, times, logfreqs] = std_ersp(EEG, ...
%                       'type', 'ersp', 'freqs', [3 50], ...
%                                 'cycles', [3 0.5], 'alpha', 0.01);
%
% See also: timef(), std_itc(), std_erp(), std_spec(), std_topo(), std_preclust()
%
% Authors: Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, January, 2005-

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

function [X, times, freqs, parameters] = std_ersp(EEG, varargin)

if nargin < 1
    help std_ersp;
    return;
end;

options = {};
if length(varargin) > 1 
    if ~isstr(varargin{1})
        if length(varargin) > 0, options = { options{:} 'components' varargin{1} }; end;
        if length(varargin) > 1, options = { options{:} 'freqs'      varargin{2} }; end;
        if length(varargin) > 2, options = { options{:} 'timewindow' varargin{3} }; end;
        if length(varargin) > 3, options = { options{:} 'cycles'     varargin{4} }; end;
        if length(varargin) > 4, options = { options{:} 'padratio'   varargin{5} }; end;
        if length(varargin) > 5, options = { options{:} 'alpha'      varargin{6} }; end;
        if length(varargin) > 6, options = { options{:} 'type'       varargin{7} }; end;
        if length(varargin) > 7, options = { options{:} 'powbase'    varargin{8} }; end;
    else
        options = varargin;
    end;
end;

[g timefargs] = finputcheck(options, { ...
                        'components'    'integer'     []      [];
                        'channels'      { 'cell' 'integer' }  { [] [] }     {};
                        'outputfile'    'string'      []      '';
                        'powbase'       'real'        []      [];
                        'trialindices' { 'integer' 'cell' } []         [];
                        'savetrials'    'string'      { 'on' 'off' }      'off';
                        'plot'          'string'      { 'on' 'off' }      'off';
                        'recompute'     'string'      { 'on' 'off' }      'off';
                        'getparams'     'string'      { 'on' 'off' }      'off';
                        'timewindow'    'real'        []      [];
                        'fileout'       'string'  []         '';
                        'timelimits'    'real'        []      [EEG.xmin EEG.xmax]*1000;
                        'cycles'        'real'        []      [3 .5];
                        'padratio'      'real'        []      1;
                        'freqs'         'real'        []      [0 EEG.srate/2];
                        'rmcomps'       'cell'        []      cell(1,length(EEG));
                        'interp'        'struct'      { }     struct([]);
                        'freqscale'     'string'      []      'log';
                        'alpha'         'real'        []      NaN;
                        'type'          'string'      { 'ersp' 'itc' 'both' 'ersp&itc' }  'both'}, 'std_ersp', 'ignore');
if isstr(g), error(g); end;
if isempty(g.trialindices), g.trialindices = cell(length(EEG)); end;
if ~iscell(g.trialindices), g.trialindices = { g.trialindices }; end;
    
% checking input parameters
% -------------------------
if isempty(g.components) & isempty(g.channels)
    if isempty(EEG(1).icaweights)
        error('EEG.icaweights not found');
    end
    g.components = 1:size(EEG(1).icaweights,1);
    disp('Computing ERSP with default values for all components of the dataset');
end

% select ICA components or data channels
% --------------------------------------
if isempty(g.fileout), g.fileout = fullfile(EEG(1).filepath, EEG(1).filename(1:end-4)); end;
if ~isempty(g.components)
    g.indices = g.components;
    prefix = 'comp';
    filenameersp   = [ g.fileout '.icaersp'  ];
    filenameitc    = [ g.fileout '.icaitc'   ];
    filenametrials = [ g.fileout '.icatimef' ];    
    if ~isempty(g.channels)
        error('Cannot compute ERSP/ITC for components and channels at the same time');
    end;
elseif ~isempty(g.channels)
    g.indices = g.channels;
    prefix = 'chan';
    filenameersp   = [ g.fileout '.datersp'  ];
    filenameitc    = [ g.fileout '.datitc'   ];
    filenametrials = [ g.fileout '.dattimef' ];    
end;

powbaseexist = 1; % used also later
if isempty(g.powbase) | isnan(g.powbase)
    powbaseexist = 0;
    g.powbase = NaN*ones(length(g.indices),1);  % default for timef()
end;
if size(g.powbase,1) ~= length(g.indices)
    error('powbase should be of size (ncomps,nfreqs)');
end

% Check if ERSP/ITC information found in datasets and if fits requested parameters 
% ----------------------------------------------------------------------------
if exist( filenameersp ) & strcmpi(g.recompute, 'off')
    fprintf('Use existing file for ERSP: %s\n', filenameersp);
    return;
end;
%    tmpersp  = load( '-mat', filenameersp, 'parameters'); % AND IT SHOULD BE USED HERE TOO - ARNO
%	params   = struct(tmpersp.parameters{:});
%    if ~isequal(params.cycles, g.cycles)                   ...
%            | (g.padratio ~= params.padratio) ...
%            | ( (g.alpha~= params.alpha) & ~( isnan(g.alpha) & isnan(params.alpha)) )
%        % if not computed with the requested parameters, recompute ERSP/ITC
%        % i.e., continue
%    else
%        disp('File ERSP/ITC data already present, computed with the same parameters: no need to recompute...');
%        return; % no need to compute ERSP/ITC
%    end
%end;

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
end;
if powbaseexist & time_range(1) >= 0 
    parameters{end+1} = 'baseboot';
    parameters{end+1} = 0;
    fprintf('No pre-0 baseline spectral estimates: Using whole epoch for timef() "baseboot"\n');
end

% return parameters
% -----------------
if strcmpi(g.getparams, 'on')
    X = []; times = []; freqs = [];
    if strcmpi(g.savetrials, 'on')
        parameters = { parameters{:} 'savetrials', g.savetrials };
    end;
    return;
end;

% No usable ERSP/ITC information available
% ---------------------------------
% tmpdata = [];
% for index = 1:length(EEG)
%     if isstr(EEG(index).data)
%         TMP = eeg_checkset( EEG(index), 'loaddata' );  % load EEG.data and EEG.icaact
%     else
%         TMP = EEG;
%     end
%     if ~isempty(g.components)
%         if isempty(TMP.icaact)                      % make icaact if necessary
%             TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
%                           reshape(TMP.data(TMP.icachansind,:,:), [ length(TMP.icachansind) size(TMP.data,2)*size(TMP.data,3) ]);
%         end;
%         tmpdata    = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2) size(TMP.data,3) ]);
%         tmpdata    = tmpdata(g.indices, :,:);
%     else
%         if isempty(tmpdata)
%             tmpdata = TMP.data(g.indices,:,:);
%         else    
%             tmpdata(:,:,end+1:end+size(TMP.data,3)) = TMP.data(g.indices,:,:);
%         end;
%     end;
% end;

options = {};
if ~isempty(g.rmcomps), options = { options{:} 'rmcomps' g.rmcomps }; end;
if ~isempty(g.interp),  options = { options{:} 'interp' g.interp }; end;
if ~isempty(g.components)
    tmpdata = eeg_getdatact(EEG, 'component', g.indices, 'trialindices', g.trialindices{1});
else
    EEG.data = eeg_getdatact(EEG, 'channel', [1:EEG.nbchan], 'rmcomps', g.rmcomps{1}, 'trialindices', g.trialindices{1});
    EEG.trials = size(EEG.data,3);
    EEG.event  = [];
    EEG.epoch  = [];
    
    if ~isempty(g.interp), 
        EEG = eeg_interp(EEG, g.interp, 'spherical'); 
    end;
    
    % find channel index
    % ------------------
    if iscell(g.channels)
        for index = 1:length(g.channels)
            chanind = strmatch( lower(g.channels{index}), lower({ EEG.chanlocs.labels }), 'exact');
            if isempty(chanind), error('Channel group not found'); end;
            chaninds(index) = chanind;
        end;
        g.indices  = chaninds;
        g.channels = chaninds;
    end;
    
    if ~iscell(g.indices)
        tmpdata = EEG.data(g.indices,:,:);
    end;
end;

% frame range
% -----------
pointrange1 = round(max((g.timelimits(1)/1000-EEG(1).xmin)*EEG(1).srate, 1));
pointrange2 = round(min(((g.timelimits(2)+1000/EEG.srate)/1000-EEG(1).xmin)*EEG(1).srate, EEG(1).pnts));
pointrange = [pointrange1:pointrange2];

% Compute ERSP & ITC
% ------------------
all_ersp   = [];
all_trials = [];
all_itc    = [];
for k = 1:length(g.indices)  % for each (specified) component
    if k>size(tmpdata,1), break; end; % happens for components
    if powbaseexist
        tmpparams = parameters;
        tmpparams{end+1} = 'powbase';
        tmpparams{end+1} = g.powbase(k,:);
    else
        tmpparams = parameters;
    end;
    
    % Run timef() to get ERSP
    % ------------------------
    timefdata  = reshape(tmpdata(k,pointrange,:), 1, length(pointrange)*size(tmpdata,3));
    if strcmpi(g.plot, 'on'), figure; end;
    [logersp,logitc,logbase,times,logfreqs,logeboot,logiboot,alltfX] ...
          = newtimef( timefdata, length(pointrange), g.timelimits, EEG(1).srate, tmpparams{2:end});
    %figure; newtimef( TMP.data(32,:), EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles, 'freqs', freqs);
    %figure; newtimef( timefdata, length(pointrange), g.timelimits, EEG.srate, cycles, 'freqs', freqs);
    if strcmpi(g.plot, 'on'), return; end;

    all_ersp = setfield( all_ersp, [ prefix int2str(g.indices(k)) '_ersp'     ], single(logersp ));
    all_ersp = setfield( all_ersp, [ prefix int2str(g.indices(k)) '_erspbase' ], single(logbase ));
    all_ersp = setfield( all_ersp, [ prefix int2str(g.indices(k)) '_erspboot' ], single(logeboot));
    all_itc  = setfield( all_itc , [ prefix int2str(g.indices(k)) '_itc'      ], single(logitc  ));
    all_itc  = setfield( all_itc , [ prefix int2str(g.indices(k)) '_itcboot'  ], single(logiboot));

    if strcmpi(g.savetrials, 'on')
        all_trials = setfield( all_trials, [ prefix int2str(g.indices(k)) '_timef'     ], single( alltfX ));
    end;
end

% Save ERSP into file
% -------------------
all_ersp.freqs      = logfreqs;
all_ersp.times      = times;
all_ersp.datatype   = 'ERSP';
all_itc.freqs       = logfreqs;
all_itc.times       = times;
all_itc.parameters  = parameters;
all_itc.datatype    = 'ITC';
all_trials.freqs    = logfreqs;
all_trials.times    = times;
all_trials.parameters = { options{:} parameters{:} };
all_trials.datatype   = 'TIMEF';

if powbaseexist
    all_ersp.parameters = { parameters{:}, 'baseline', g.powbase };
else
    all_ersp.parameters = parameters;
end;
if ~isempty(g.channels)
    if ~isempty(EEG(1).chanlocs)
        all_ersp.chanlabels   = { EEG(1).chanlocs(g.channels).labels };
        all_itc.chanlabels    = { EEG(1).chanlocs(g.channels).labels };
        all_trials.chanlabels = { EEG(1).chanlocs(g.channels).labels };
    end;
end;

if strcmpi(g.type, 'both') | strcmpi(g.type, 'ersp') | strcmpi(g.type, 'ersp&itc')
    std_savedat( filenameersp, all_ersp);
end;
if strcmpi(g.type, 'both') | strcmpi(g.type, 'itc') | strcmpi(g.type, 'ersp&itc')
    std_savedat( filenameitc , all_itc );
end;
if strcmpi(g.savetrials, 'on')
    std_savedat( filenametrials , all_trials );
end;

