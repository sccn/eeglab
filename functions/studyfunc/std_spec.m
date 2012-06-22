% std_spec() - Returns the ICA component spectra for a dataset. Updates the EEG structure 
%              in the Matlab environment and in the .set file as well. Saves the spectra 
%              in a file.
% Usage:    
%           >> [spec freqs] = std_spec(EEG, 'key', 'val', ...);
%
%              Computes the mean spectra of the activites of specified components of the 
%              supplied dataset. The spectra are saved in a Matlab file. If such a file 
%              already exists, loads the spectral information from this file.  
%              Options (below) specify which components to use, and the desired frequency 
%              range. There is also an option to specify other spectopo() input variables 
%              (see >> help spectopo for details).
%
%              Returns the removed mean spectra of the selected ICA components in the 
%              requested frequency range. If the spectra were computed previously but a
%              different frequency range is selected, there is an overwrite option. 
%              so. The function will load previously computed log spectra, if any, and 
%              will remove the mean from the requested frequency range. The frequencies 
%              vector is also returned. 
% Inputs:
%   EEG - a loaded epoched EEG dataset structure. 
%
% Optional inputs:
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
%   'trialindices' - [cell array] indices of trials for each dataset.
%                  Default is all trials.
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
%                  single-trial statistics.
%
% spectrum specific optional inputs:
%   'specmode'   - ['psd'|'fft'|'pburg'|'pmtm'] method to compute spectral 
%                  decomposition. 'psd' uses the spectopo function. 'fft' 
%                  uses a simple fft on each trial. For continuous data
%                  data trials are extracted automatically (see 'epochlim'
%                  and 'epochrecur' below). Two experimental modes are 
%                  'pmtm' and 'pbug' which use multitaper and the Burg 
%                  method to compute spectrum respectively.
%   'epochlim'   - [min max] for FFT on continuous data, extract data
%                  epochs with specific epoch limits in seconds (see also
%                  'epochrecur' below). Default is [0 1].
%   'epochrecur' - [float] for FFT on continuous data, set the automatic
%                  epoch extraction recurence interval (default is 0.5 second).
%   'timerange'  - [min max] use data within a specific time range before 
%                  computing the data spectrum. For instance, for evoked 
%                  data trials, it is recommended to use the baseline time 
%                  period.
%   'continuous' - ['on'|'off'] force epoch data to be treated as
%                  continuous so small data epochs can be extracted for the
%                  'fft' specmode option. Default is 'off'.
%   'freqrange'  - [minhz maxhz] frequency range (in Hz) within which to 
%                  return the spectrum {default|[]: [0 sample rate/2]}.
%                  Note that this does not affect the spectrum computed on
%                  disk, only the data returned by this function as output.
%   'nw'         - [integer] number of tapers for the 'pmtm' spectral
%                  method. Default is 4.
%   'burgorder'  - [integet] order for the Burg spectral method.
%
% Other optional spectral parameters:
%   All optional parameters to the spectopo function may be provided to this 
%   function as well (requires the 'specmode' option above to be set to
%   'psd').
%
% Outputs:
%   spec      - the mean spectra (in dB) of the requested ICA components in the selected 
%               frequency range (with the mean of each spectrum removed). 
%   freqs     - a vector of frequencies at which the spectra have been computed. 
%
% Files output or overwritten for ICA: 
%               [dataset_filename].icaspec,   % raw spectrum of ICA components
% Files output or overwritten for data: 
%               [dataset_filename].datspec, 
% 
% See also  spectopo(), std_erp(), std_ersp(), std_map(), std_preclust()
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, January, 2005

% Defunct:      0 -> if frequency range is different from saved spectra, ask via a 
%                    pop-up window whether to keep existing spectra or to overwrite them. 

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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

function [X, f, overwrt] = std_spec(EEG, varargin)

overwrt = 1; % deprecated
if nargin < 1
    help std_spec;
    return;
end;

% decode inputs
% -------------
if ~isempty(varargin) 
    if ~isstr(varargin{1})
        varargin = { varargin{:} [] [] };
        if all(varargin{1} > 0) 
            options = { 'components' varargin{1} 'freqrange' varargin{2} };
        else
            options = { 'channels' -varargin{1} 'freqrange' varargin{2} };
        end;
    else
        options = varargin;
    end;
else
    options = varargin;
end;

[g spec_opt] = finputcheck(options, { 'components' 'integer' []         [];
                                      'channels'   'cell'    {}         {};
                                      'timerange'  'float'   []         [];
                                      'specmode'   'string'  {'fft','psd','pmtm','pburg'} 'psd';
                                      'recompute'  'string'  { 'on','off' } 'off';
                                      'savetrials' 'string'  { 'on','off' } 'off';
                                      'continuous' 'string'  { 'on','off' } 'off';
                                      'savefile'   'string'  { 'on','off' } 'on';
                                      'epochlim'   'real'    []         [0 1];
                                      'trialindices' { 'integer','cell' } []         [];
                                      'epochrecur' 'real'    []         0.5;
                                      'rmcomps'    'cell'    []         cell(1,length(EEG));
                                      'nw'         'float'   []         4;
                                      'fileout'    'string'  []         '';
                                      'burgorder'  'integer' []         20;
                                      'interp'     'struct'  { }        struct([]);
                                      'nfft'       'integer' []         [];
                                      'freqrange'  'real'    []         [] }, 'std_spec', 'ignore');
if isstr(g), error(g); end;
if isfield(EEG,'icaweights')
   numc = size(EEG(1).icaweights,1);
else
   error('EEG.icaweights not found');
end
if isempty(g.components)
    g.components = 1:numc;
end

EEG_etc = [];

% filename 
% --------
if isempty(g.fileout), g.fileout = fullfile(EEG(1).filepath, EEG(1).filename(1:end-4)); end;
if ~isempty(g.channels)
    filename = [ g.fileout '.datspec'];
    prefix = 'chan';
else    
    filename = [ g.fileout '.icaspec'];
    prefix = 'comp';
end;

% SPEC information found in datasets
% ---------------------------------
if exist(filename) & strcmpi(g.recompute, 'off')

    fprintf('File "%s" found on disk, no need to recompute\n', filename);
    setinfo.filebase = g.fileout;
    if strcmpi(prefix, 'comp')
        [X tmp f] = std_readfile(setinfo, 'components', g.components, 'freqlimits', g.freqrange, 'measure', 'spec');
    else
        [X tmp f] = std_readfile(setinfo, 'channels', g.channels, 'freqlimits', g.freqrange, 'measure', 'spec');
    end;
    if ~isempty(X), return; end;
end

% No SPEC information found
% -------------------------
options = {};
if ~isempty(g.rmcomps), options = { options{:} 'rmcomps' g.rmcomps }; end;
if ~isempty(g.interp),  options = { options{:} 'interp' g.interp }; end;
if isempty(g.channels)
     [X boundaries]  = eeg_getdatact(EEG, 'component', [1:size(EEG(1).icaweights,1)], 'trialindices', g.trialindices );
else [X boundaries]  = eeg_getdatact(EEG, 'channel'  , [1:EEG(1).nbchan], 'trialindices', g.trialindices, 'rmcomps', g.rmcomps, 'interp', g.interp);
end;
if ~isempty(boundaries) && boundaries(end) ~= size(X,2), boundaries = [boundaries size(X,2)]; end;

% extract epochs if necessary
% ---------------------------
oritrials = EEG.trials;
if ~strcmpi(g.specmode, 'psd')
    if EEG(1).trials == 1 || strcmpi(g.continuous, 'on')
        EEG.data = X;
        EEG.trials = size(EEG.data,3);
        if EEG(1).trials > 1
            EEG = eeg_epoch2continuous(EEG);
        end;
        EEG = eeg_regepochs(EEG, g.epochrecur, g.epochlim);
        g.trialindices = { [1:EEG(1).trials] };
        disp('Warning: continuous data, extracting 1-second epochs'); 
        X = EEG.data;
    end;
end;
 
% get specific time range for epoched and continuous data
if ~isempty(g.timerange) 
    if oritrials > 1
        timebef  = find(EEG(1).times >= g.timerange(1) & EEG(1).times < g.timerange(2) );
        X        = X(:,timebef,:);
        EEG(1).pnts = length(timebef);
    else
        disp('warning: ''timerange'' option cannot be used with continuous data');
    end;
end;

% compute spectral decomposition
if strcmpi(g.specmode, 'psd')
    [X, f] = spectopo(X, size(X,2), EEG(1).srate, 'plot', 'off', 'boundaries', boundaries, 'nfft', g.nfft, spec_opt{:});
    if strcmpi(g.savetrials, 'on')
        disp('WARNING: Cannot save trials using ''psd'' specmode option');
    end;
elseif strcmpi(g.specmode, 'pmtm')
    if all([ EEG.trials ] == 1) && ~isempty(boundaries), disp('Warning: multitaper does not take into account boundaries in continuous data'); end;
    fprintf('Computing multitaper:');
    for cind = 1:size(X,1)
        fprintf('.');
        for tind = 1:size(X,3)
            [tmpdat f] = pmtm(X(cind,:,tind), g.nw, g.nfft, EEG.srate);
            if cind == 1 && tind == 1
                X2 = zeros(size(X,1), length(tmpdat), size(X,3));
            end;
            X2(cind,:,tind) = tmpdat;
        end;
    end;
    fprintf('\n');
    X     = 10*log10(X2);  
    if strcmpi(g.savetrials, 'off')
        X = mean(X,3);
    end;
elseif strcmpi(g.specmode, 'pburg')
    if all([ EEG.trials ] == 1) && ~isempty(boundaries), disp('Warning: pburg does not take into account boundaries in continuous data'); end;
    for cind = 1:size(X,1)
        fprintf('.');
        for tind = 1:size(X,3)
            [tmpdat f] = pburg(X(cind,:,tind), g.burgorder, g.nfft, EEG.srate);
            if cind == 1 && tind == 1
                X2 = zeros(size(X,1), length(tmpdat), size(X,3));
            end;
            X2(cind,:,tind) = tmpdat;
        end;
    end;
    fprintf('\n');
    if strcmpi(g.savetrials, 'off')
        X = mean(X,3);    
    end;
else % fft mode
    if oritrials == 1 || strcmpi(g.continuous, 'on')
        X = bsxfun(@times, X, hamming(size(X,2))');
    end;
    if all([ EEG.trials ] == 1) && ~isempty(boundaries), disp('Warning: fft does not take into account boundaries in continuous data'); end;
    tmp   = fft(X, g.nfft, 2);
    f     = linspace(0, EEG(1).srate/2, floor(size(tmp,2)/2));
    f     = f(2:end); % remove DC (match the output of PSD)
    tmp   = tmp(:,2:floor(size(tmp,2)/2),:);
    X     = 10*log10(abs(tmp).^2);    
    if strcmpi(g.savetrials, 'off')
        X = mean(X,3);  
    end;
end;

% Save SPECs in file (all components or channels)
% ----------------------------------
fileNames = computeFullFileName( { EEG.filepath }, { EEG.filename });
if strcmpi(g.savefile, 'on')
    options = { options{:} spec_opt{:} 'timerange' g.timerange 'nfft' g.nfft 'specmode' g.specmode };
    if strcmpi(prefix, 'comp')
        savetofile( filename, f, X, 'comp', 1:size(X,1), options, {}, fileNames, g.trialindices);
    else
        if ~isempty(g.interp)
            savetofile( filename, f, X, 'chan', 1:size(X,1), options, { g.interp.labels }, fileNames, g.trialindices);
        else
            tmpchanlocs = EEG(1).chanlocs;
            savetofile( filename, f, X, 'chan', 1:size(X,1), options, { tmpchanlocs.labels }, fileNames, g.trialindices);
        end;
    end;
end;
return;

% compute full file names
% -----------------------
function res = computeFullFileName(filePaths, fileNames);
for index = 1:length(fileNames)
    res{index} = fullfile(filePaths{index}, fileNames{index});
end;

% -------------------------------------
% saving SPEC information to Matlab file
% -------------------------------------
function savetofile(filename, f, X, prefix, comps, params, labels, dataFiles, dataTrials);
    
    disp([ 'Saving SPECTRAL file ''' filename '''' ]);
    allspec = [];
    for k = 1:length(comps)
        allspec = setfield( allspec, [ prefix int2str(comps(k)) ], squeeze(X(k,:,:)));
    end;
    if nargin > 6 && ~isempty(labels)
        allspec.labels = labels;
    end;
    allspec.freqs      = f;
    allspec.parameters = params;
    allspec.datatype   = 'SPECTRUM';
    allspec.datafiles   = dataFiles;
    allspec.datatrials  = dataTrials;
    allspec.average_spec = mean(X,1);
    std_savedat(filename, allspec);

