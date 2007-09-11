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

% $Log: not supported by cvs2svn $
% Revision 1.53  2007/08/12 23:20:24  arno
% clarifying help message
%
% Revision 1.52  2007/08/12 04:10:33  arno
% help message plus allow g.channels to be an integer array
%
% Revision 1.51  2007/08/09 22:12:57  arno
% file exist for ERSP
%
% Revision 1.50  2007/05/22 13:56:32  arno
% default highest frequency
%
% Revision 1.49  2007/04/27 19:05:45  julie
% added tmpparams to use input baseline
%
% Revision 1.48  2007/04/27 18:49:08  arno
% baseline
%
% Revision 1.47  2007/04/07 21:40:33  arno
% ersp and itc
%
% Revision 1.46  2007/04/05 23:13:43  arno
% recomputing file on disk
%
% Revision 1.45  2007/04/05 22:17:02  arno
% freqscale
%
% Revision 1.44  2007/04/05 22:09:38  arno
% default for freqscale
%
% Revision 1.43  2007/04/05 20:59:46  arno
% default is log for frequencies
%
% Revision 1.42  2006/11/14 03:11:24  arno
% nothing
%
% Revision 1.41  2006/10/05 15:13:28  arno
% fixing processing ICA data for 'savetrials'
%
% Revision 1.40  2006/10/03 13:38:22  scott
% worked on help msg. ARNO - See many ?? in the help text! -sm
%
% Revision 1.39  2006/10/02 11:40:51  arno
% minor changes
%
% Revision 1.35  2006/04/11 18:24:07  arno
% fixing alpha NaN problem
%
% Revision 1.34  2006/03/16 02:52:27  scott
% adding 'baseline',powbase to saved params -sm
%
% Revision 1.33  2006/03/16 01:46:12  scott
% fix appending to cell array parameters
%
% Revision 1.32  2006/03/16 01:12:49  scott
% added back? comma
%
% Revision 1.31  2006/03/15 22:52:44  scott
% if powbase specified and no pre-0 spectral estimates, use whole epoch for baseboot (with warning)
%
% Revision 1.30  2006/03/15 19:49:07  scott
% nothing
%
% Revision 1.29  2006/03/12 02:51:02  arno
% function
% call
%
% Revision 1.28  2006/03/11 07:11:28  arno
% header
%
% Revision 1.27  2006/03/11 00:35:50  arno
% retreive previous version
%std_filecheck
% Revision 1.25  2006/03/10 15:50:07  arno
% converting values to single
%
% Revision 1.24  2006/03/10 00:30:55  arno
% update header
%
% Revision 1.23  2006/03/09 23:29:21  arno
% implement new ERSP from Matlab and different structure ec...
%
% Revision 1.22  2006/03/09 00:37:55  arno
% nothing yet
%
% Revision 1.21  2006/03/08 23:03:25  arno
% remove debug message
%
% Revision 1.20  2006/03/08 23:02:41  arno
% renaming powbase properly
%
% Revision 1.19  2006/03/08 22:55:50  arno
% fix itcboot now
%
% Revision 1.18  2006/03/08 22:50:40  arno
% fix last change
%
% Revision 1.17  2006/03/08 22:49:10  arno
% detect the presence of file or not
%
% Revision 1.16  2006/03/08 20:29:08  arno
% rename func
%
% Revision 1.15  2006/03/08 19:43:37  scott
% fixed powbase definition -sm & ad
%
% Revision 1.14  2006/03/08 03:02:27  scott
% expand powbase to a matrix -sm
%
% Revision 1.13  2006/03/08 02:58:08  scott
% debug3
%
% Revision 1.12  2006/03/08 02:54:41  scott
% debug 2
%
% Revision 1.11  2006/03/08 02:38:14  scott
% debug
%
% Revision 1.10  2006/03/07 23:03:47  scott
% added optional posbase argument for timef(). NOTE: probability masking now
% allowed and computed
% but the masked version is not output by timef (I think - must clarify timef help)
% and therefore not used here... -sm
%
% Revision 1.9  2006/03/07 19:18:44  arno
% header
%
% Revision 1.8  2006/03/07 03:58:32  scott
% same -sm
%
% Revision 1.7  2006/03/07 03:32:22  scott
% fixing default|specified comps calculation -sm
%
% Revision 1.6  2006/03/07 02:30:54  scott
% worked on help msg; formatted, made accurate (only 2 files are now output, not 4).
% made the function accept the components argument (was always computing ersp/itc
% for ALL components, not just those asked for.  -sm
%
% Revision 1.5  2006/03/06 23:18:08  arno
% resave
%
% Revision 1.4  2006/03/03 23:39:05  arno
% put log
%

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
                        'savetrials'    'string'      { 'on' 'off' }      'off';
                        'plot'          'string'      { 'on' 'off' }      'off';
                        'recompute'     'string'      { 'on' 'off' }      'off';
                        'getparams'     'string'      { 'on' 'off' }      'off';
                        'timewindow'    'real'        []      [];
                        'timelimits'    'real'        []      [EEG.xmin EEG.xmax]*1000;
                        'cycles'        'real'        []      [3 .5];
                        'padratio'      'real'        []      1;
                        'freqs'         'real'        []      [3 EEG.srate/2];
                        'freqscale'     'string'      []      'log';
                        'alpha'         'real'        []      NaN;
                        'type'          'string'      { 'ersp' 'itc' 'both' 'ersp&itc' }  'both'}, 'std_ersp', 'ignore');
if isstr(g), error(g); end;
    
% checking input parameters
% -------------------------
if isempty(g.components) & isempty(g.channels)
    if isempty(EEG(1).icaweights)
        error('EEG.icaweights not found');
    end
    g.components = 1:size(EEG(1).icaweights,1);
    disp('Computing ERSP with default values for all components of the dataset');
end

% find channel index
% ------------------
if ~isempty(g.channels)
    if iscell(g.channels)
        for index = 1:length(g.channels)
            chanind = strmatch( lower(g.channels{index}), lower({ EEG.chanlocs.labels }), 'exact');
            if isempty(chanind), error('Channel group not found'); end;
            chaninds(index) = chanind;
        end;
        g.channels = chaninds;
    end;
end;

% select ICA components or data channels
% --------------------------------------
if ~isempty(g.outputfile)
    filenameersp   = fullfile('', [ g.outputfile '.datersp' ]);
    filenameitc    = fullfile('', [ g.outputfile '.datitc' ]);    
    filenametrials = fullfile('', [ g.outputfile '.dattimef' ]);    
    g.indices = g.channels;
    prefix = 'chan';
elseif ~isempty(g.components)
    g.indices = g.components;
    prefix = 'comp';
    filenameersp   = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaersp' ]);
    filenameitc    = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaitc' ]);
    filenametrials = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icatimef' ]);    
    if ~isempty(g.channels)
        error('Cannot compute ERSP/ITC for components and channels at the same time');
    end;
elseif ~isempty(g.channels)
    g.indices = g.channels;
    prefix = 'chan';
    filenameersp   = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'datersp' ]);
    filenameitc    = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'datitc' ]);
    filenametrials = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'dattimef' ]);    
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
[time_range, g.winsize] = compute_ersp_times(g.cycles,  EEG(1).srate, ...
                                 [EEG(1).xmin EEG(1).xmax]*1000 , g.freqs(1), g.padratio); 
if time_range(1) >= time_range(2)
    error(['std_ersp(): parameters given for ' upper(g.type) ...
           'calculation result in an invalid time range. Aborting.' ...
           'Please increase the lower frequency bound or change other' ...
           'parameters to resolve the problem. See >> timef details'] )
end
parameters = { 'cycles', g.cycles, 'padratio', g.padratio,  'winsize', round(g.winsize), ...
               'alpha', g.alpha, 'freqscale', g.freqscale, timefargs{:} };
if length(g.freqs)>0, parameters = { parameters{:} 'freqs' g.freqs }; end;
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
    return;
end;

% No usable ERSP/ITC information available
% ---------------------------------
tmpdata = [];
for index = 1:length(EEG)
    if isstr(EEG(index).data)
        TMP = eeg_checkset( EEG(index), 'loaddata' );  % load EEG.data and EEG.icaact
    else
        TMP = EEG;
    end
    if ~isempty(g.components)
        if isempty(TMP.icaact)                      % make icaact if necessary
            TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
                          reshape(TMP.data(TMP.icachansind,:,:), [ length(TMP.icachansind) size(TMP.data,2)*size(TMP.data,3) ]);
        end;
        tmpdata    = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2) size(TMP.data,3) ]);
        tmpdata    = tmpdata(g.indices, :,:);
    else
        if isempty(tmpdata)
            tmpdata = TMP.data(g.indices,:,:);
        else    
            tmpdata(:,:,end+1:end+size(TMP.data,3)) = TMP.data(g.indices,:,:);
        end;
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
all_trials.parameters = parameters;
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

