% std_erp() -   Constructs and returns channel or ICA activation ERPs for a dataset. 
%               Saves the ERPs into a Matlab file, [dataset_name].icaerp, for
%               data channels or [dataset_name].icaerp for ICA components, 
%               in the same directory as the dataset file.  If such a file 
%               already exists, loads its information. 
% Usage:    
%            >> [erp, times] = std_erp(EEG, 'key', 'val', ...);
% Inputs:
%   EEG          - a loaded epoched EEG dataset structure. May be an array
%                  of such structure containing several datasets.
%
% Optional inputs:
%   'components' - [numeric vector] components of the EEG structure for which 
%                  activation ERPs will be computed. Note that because 
%                  computation of ERP is so fast, all components ERP are
%                  computed and saved. Only selected component 
%                  are returned by the function to Matlab
%                  {default|[] -> all}
%   'channels'   - [cell array] channels of the EEG structure for which 
%                  activation ERPs will be computed. Note that because 
%                  computation of ERP is so fast, all channels ERP are
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
%                  to interpolate (this entry is ignored when plotting 
%                  components). Default is no interpolation.
%   'fileout'    - [string] name of the file to save on disk. The default
%                  is the same name (with a different extension) as the 
%                  dataset given as input.
%  'savetrials'  - ['on'|'off'] save single-trials ERSP. Requires a lot of disk
%                  space (dataset space on disk times 10) but allow for refined
%                  single-trial statistics.
%
% ERP specific options:
%   'rmbase'     - [min max] remove baseline. This option does not affect
%                  the original datasets.
%
% Outputs:
%   erp          - ERP for the requested ICA components in the selected 
%                  latency window. ERPs are scaled by the RMS over of the
%                  component scalp map projection over all data channels.
%   times        - vector of times (epoch latencies in ms) for the ERP
%
% File output:     
%    [dataset_file].icaerp     % component erp file
% OR
%    [dataset_file].daterp     % channel erp file
%
% See also: std_spec(), std_ersp(), std_topo(), std_preclust()
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, January, 2005

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

function [X, t] = std_erp(EEG, varargin); %comps, timerange)

if nargin < 1
    help std_erp;
    return;
end

% decode inputs
% -------------
if ~isempty(varargin) 
    if ~ischar(varargin{1})
        varargin = { varargin{:} [] [] };
        if all(varargin{1} > 0) 
            options = { 'components' varargin{1} 'timerange' varargin{2} };
        else
            options = { 'channels' -varargin{1} 'timerange' varargin{2} };
        end
    else
        options = varargin;
    end
else
    options = varargin;
end

g = finputcheck(options, { 'components' 'integer' []         [];
                           'channels'   'cell'    {}         {};
                           'rmbase'     'real'    []         [];
                           'trialindices' { 'integer','cell' } []         [];
                           'rmcomps'    'cell'    []         cell(1,length(EEG));
                           'fileout'    'string'  []         '';
                           'trialinfo'  'struct'  []         struct([]);
                           'savetrials' 'string'  { 'on','off' } 'off';
                           'interp'     'struct'  { }        struct([]);
                           'timerange'  'real'    []         [];        % the timerange option is deprecated and has no effect
                           'recompute'  'string'  { 'on','off' } 'off' }, 'std_erp');
if ischar(g), error(g); end
if isempty(g.trialindices), g.trialindices = cell(length(EEG)); end
if ~iscell(g.trialindices), g.trialindices = { g.trialindices }; end
if isfield(EEG,'icaweights')
   numc = size(EEG(1).icaweights,1);
else
   error('EEG.icaweights not found');
end
if isempty(g.components)
    g.components = 1:numc;
end

% % THIS SECTION WOULD NEED TO TEST THAT THE PARAMETERS ON DISK ARE CONSISTENT
%
% % filename 
% % --------
if isempty(g.fileout), g.fileout = fullfile(EEG(1).filepath, EEG(1).filename(1:end-4)); end
if ~isempty(g.channels)
    filenameshort = [ g.fileout '.daterp'];
    prefix = 'chan';
else    
    filenameshort = [ g.fileout '.icaerp'];
    prefix = 'comp';
end
%filename = fullfile( EEG(1).filepath, filenameshort);
filename = filenameshort;

% ERP information found in datasets
% ---------------------------------
if exist(filename) && strcmpi(g.recompute, 'off') && nargout > 0

    fprintf('File "%s" found on disk, no need to recompute\n', filenameshort);
    setinfo.filebase = filename;
    if strcmpi(prefix, 'comp')
        [X tmp t] = std_readfile(setinfo, 'components', g.components, 'timelimits', g.timerange, 'measure', 'erp');
    else
        [X tmp t] = std_readfile(setinfo, 'channels', g.channels,  'timelimits', g.timerange, 'measure', 'erp');
    end
    if ~isempty(X), return; end
    
end 
   
% No ERP information found
% ------------------------
% if ischar(EEG.data)
%     TMP = eeg_checkset( EEG, 'loaddata' ); % load EEG.data and EEG.icaact
% else
%     TMP = EEG;
% end
%    && isempty(TMP.icaact)
%    TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
%        reshape(TMP.data(TMP.icachansind,:,:), [ length(TMP.icachansind) size(TMP.data,2)*size(TMP.data,3) ]);
%    TMP.icaact = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2) size(TMP.data,3) ]);
%end
%if strcmpi(prefix, 'comp'), X = TMP.icaact;
%else                        X = TMP.data;
%end
if isempty(g.channels)
     X = eeg_getdatact(EEG, 'component', [1:size(EEG(1).icaweights,1)], 'trialindices', g.trialindices );
else X = eeg_getdatact(EEG, 'trialindices', g.trialindices, 'rmcomps', g.rmcomps, 'interp', g.interp);
end

% Remove baseline mean
% --------------------
pnts     = EEG(1).pnts;
trials   = size(X,3);
timevals = EEG(1).times;
if ~isempty(g.timerange)
    disp('Warning: the ''timerange'' option is deprecated and has no effect');
end
if ~isempty(X)
    if ~isempty(g.rmbase)
        disp('Removing baseline...');
        options = { options{:} 'rmbase' g.rmbase };
        [tmp timebeg] = min(abs(timevals - g.rmbase(1)));
        [tmp timeend] = min(abs(timevals - g.rmbase(2)));
        if ~isempty(timebeg)
            X = rmbase(X,pnts, [timebeg:timeend]);
        else
            X = rmbase(X,pnts);
        end
    end
    X = reshape(X, [ size(X,1) pnts trials ]);
    if strcmpi(prefix, 'comp')
        if strcmpi(g.savetrials, 'on')
            X = repmat(sqrt(mean(EEG(1).icawinv.^2))', [1 EEG(1).pnts size(X,3)]) .* X;
        else
            X = repmat(sqrt(mean(EEG(1).icawinv.^2))', [1 EEG(1).pnts]) .* mean(X,3); % calculate ERP
        end
    elseif strcmpi(g.savetrials, 'off')
        X = mean(X, 3);
    end
end
eeglab_options;
if option_single
    X = single(X);
end

% Save ERPs in file (all components or channels)
% ----------------------------------------------
if isempty(timevals), timevals = linspace(EEG(1).xmin, EEG(1).xmax, EEG(1).pnts)*1000; end; % continuous data
fileNames = computeFullFileName( { EEG.filepath }, { EEG.filename });
if strcmpi(prefix, 'comp')
    savetofile( filename, timevals, X, 'comp', 1:size(X,1), options, {}, fileNames, g.trialindices, g.trialinfo);
    %[X,t] = std_readerp( EEG, 1, g.components, g.timerange);
else
    if ~isempty(g.interp)
        savetofile( filename, timevals, X, 'chan', 1:size(X,1), options, { g.interp.labels }, fileNames, g.trialindices, g.trialinfo);
    else
        tmpchanlocs = EEG(1).chanlocs;
        savetofile( filename, timevals, X, 'chan', 1:size(X,1), options, { tmpchanlocs.labels }, fileNames, g.trialindices, g.trialinfo);
    end
    %[X,t] = std_readerp( EEG, 1, g.channels, g.timerange);
end

% compute full file names
% -----------------------
function res = computeFullFileName(filePaths, fileNames);
for index = 1:length(fileNames)
    res{index} = fullfile(filePaths{index}, fileNames{index});
end

% -------------------------------------
% saving ERP information to Matlab file
% -------------------------------------
function savetofile(filename, t, X, prefix, comps, params, labels, dataFiles, dataTrials, trialInfo);
    
    disp([ 'Saving ERP file ''' filename '''' ]);
    allerp = [];
    for k = 1:length(comps)
        allerp = setfield( allerp, [ prefix int2str(comps(k)) ], squeeze(X(k,:,:)));
    end
    if nargin > 6 && ~isempty(labels)
        allerp.labels = labels;
    end
    allerp.times       = t;
    allerp.datatype    = 'ERP';
    allerp.parameters  = params;
    allerp.datafiles   = dataFiles;
    allerp.datatrials  = dataTrials;
    allerp.trialinfo   = trialInfo;
    
    std_savedat(filename, allerp);
