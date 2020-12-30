% std_custom() - Returns custom measure for a dataset.
%
% Usage:    
%           >> [spec freqs] = std_custom(EEG, func, 'key', 'val', ...);
%
% Inputs:
%   EEG  - a loaded epoched EEG dataset structure
%   func - a function handle that takes a data array as input
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
%   'output'     - ['power'|'fft'] compute power of keep single complex
%                  'fft' estimate. Default is 'power'.
%   'fileout'    - [string] name of the file to save on disk. The default
%                  is the same name (with a different extension) as the 
%                  dataset given as input.
%  'savetrials'  - ['on'|'off'] save single-trials ERSP. Requires a lot of disk
%                  space (dataset space on disk times 10) but allow for refined
%                  single-trial statistics.
%
% Outputs:
%   customres - custom result array
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, Dec, 2020

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2020, arno@ucsd.edu
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

function std_custom(EEG, func, varargin)

if nargin < 1
    help std_custom;
    return;
end

% decode inputs
% -------------
options = varargin;
[g, funcopt] = finputcheck(options, { 'components' 'integer' []         [];
                                      'channels'   'cell'    {}         {};
                                      'timerange'  'float'   []         [];
                                      'recompute'  'string'  { 'on','off' } 'off';
                                      'savetrials' 'string'  { 'on','off' } 'off';
                                      'savefile'   'string'  { 'on','off' } 'on';
                                      'fileext'    'string'  { }        'custom';
                                      'trialindices' { 'integer','cell' } []         [];
                                      'rmcomps'    'cell'    []         cell(1,length(EEG));
                                      'fileout'    'string'  []         '';
                                      'trialinfo'  'struct'  []         struct([]);
                                      'fileout'    'string'  []         '';
                                      'interp'     'struct'  { }        struct([]) }, 'std_spec', 'ignore');
if ischar(g), error(g); end
if isempty(g.trialindices), g.trialindices = cell(length(EEG)); end
if ~iscell(g.trialindices), g.trialindices = { g.trialindices }; end
if isfield(EEG,'icaweights')
   numc = size(EEG(1).icaweights,1);
else
   error('EEG.icaweights not found');
end
if ~isempty(g.channels)
    filename = [ g.fileout '.dat' g.fileext];
    prefix = 'chan';
else    
    filename = [ g.fileout '.ica' g.fileext];
    prefix = 'comp';
end

if isempty(g.components)
    g.components = 1:numc;
end

% No SPEC information found
% -------------------------
if isempty(g.channels)
     [X,boundaries]  = eeg_getdatact(EEG, 'component', [1:size(EEG(1).icaweights,1)], 'trialindices', g.trialindices );
else [X,boundaries]  = eeg_getdatact(EEG, 'trialindices', g.trialindices, 'rmcomps', g.rmcomps, 'interp', g.interp);
end
if ~isempty(boundaries) && boundaries(end) ~= size(X,2), boundaries = [boundaries size(X,2)]; end
 
% get specific time range for epoched and continuous data
% -------------------------------------------------------
oritrials = EEG.trials;
if ~isempty(g.timerange) 
    if oritrials > 1
        timebef  = find(EEG(1).times >= g.timerange(1) & EEG(1).times < g.timerange(2) );
        X        = X(:,timebef,:);
        EEG(1).pnts = length(timebef);
    else
        disp('warning: ''timerange'' option cannot be used with continuous data');
    end
end

X = feval(func, X, funcopt{:});

% convert to structure
% --------------------
if strcmpi(prefix, 'comp')
    savetofile( filename, X, 'comp', 1:size(X,1), options, {}, fileNames, g.trialindices, 'custom' , g.trialinfo);
else
    if ~isempty(g.interp)
        savetofile( filename, X, 'chan', 1:size(X,1), options, { g.interp.labels }, '', g.trialindices, 'custom' , g.trialinfo);
    else
        tmpchanlocs = EEG(1).chanlocs;
        savetofile( filename, X, 'chan', 1:size(X,1), options, { tmpchanlocs.labels }, '', g.trialindices, 'custom' , g.trialinfo);
    end
end

return;

% -------------------------------------
% saving SPEC information to Matlab file
% -------------------------------------
function savetofile( filename, X, prefix, comps, params, labels, dataFiles, dataTrials, datatype , trialInfo);
    
    allcustom = [];
    for k = 1:length(comps)
        allcustom = setfield( allcustom, [ prefix int2str(comps(k)) ], squeeze(X(k,:,:)));
    end
    if nargin > 6 && ~isempty(labels)
        allcustom.labels = labels;
    end
    allcustom.parameters = params;
    allcustom.datatype   = datatype;
    allcustom.datafiles   = dataFiles;
    allcustom.datatrials  = dataTrials;
    allcustom.trialinfo   = trialInfo;
    allcustom.average_spec = mean(X,1);
    std_savedat(filename, allcustom);


