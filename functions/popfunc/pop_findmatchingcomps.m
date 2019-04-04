% pop_findmatchingcomps() - Find components with the highest 
%                absolute scalp map correlation to components provided or tagged for
%                rejection in another dataset.
%
% % Usage:
%       >> OUTEEG = pop_findmatchingcomps( INEEG, 'key', val, . . . );
%
% Inputs:
%   INEEG    - Input dataset
%
% Optional input:
%   'dataset'      - Another EEGLAB dataset. Find components in this dataset matching 
%                  those tagged for rejection in INEEG (e.g., those components
%                  marked manually in the other dataset using menu item
%                  "Tools > Reject data using ICA > Reject components by map")
%                  Else, ALLEEG index (integer) of another such dataset.
%   'corrthresh'   - [0<float<=1] minimal absolute ÔmatchingÕ correlation.
%                  Default is 0.92.
%   'matchcomps'   - [float array] scalp maps (as in icawinv) of components
%                  from another decomposition or dataset to be matched to the
%                  current EEG set decomposition.
%                  Several components may be provided as input (as different
%                  columns or different rows). Default is none.
%   'rejflag'      - [0|1] Flag to tag matching components for rejection. 
%                  Default is [0] do not tag for rejection.  
%   'nomatcherror' - ['on'|'off'] trigger error if the tagged component maps
%                  are not found.  Default is 'off'.
%
% Output:
%   OUTEEG      - Output dataset with matching artifact components tagged
%                 for rejection in OUTEEG.reject.gcompreject (if rejflag
%                 option set to 1)
%   matchic     - Indices of matching components
%   matchinput  - Indices of 'matchcomps' matching components in output
%                 mathcic. This output will be empty if option 'dataset' is
%                 used.
%
% Author: Arnaud Delorme, 2017
%
% See also: matcorr.m

% Copyright (C) 2017 Arnaud Delorme, UCSD, arno@ucsd.edu
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

% 01-25-02 reformatted help & license -ad 

function [EEG, matchic, matchinput, com] = pop_findmatchingcomps( EEG, varargin)

com = '';
if nargin < 2
    
     uilist = { { 'style'  'text'     'string'   'Index of dataset to find matching artifact comps.' } ...
               { 'style'   'edit'     'string'    '' } ...
               { 'style'   'text'     'string'    'or scalp maps (array or Matlab variable) from which to find matching artifact comps.' } ...
               { 'style'   'edit'     'string'    '' } ...
               { 'style'   'text'     'string'    '' } ...
               { 'style'   'text'     'string'    'Minimum abs. correlation for matching components' } ...
               { 'style'   'edit'     'string'    '0.92' }           };
    uigeom = { [1.5 0.5] [1] [1] [1] [1.5 0.5] };
    result = inputgui( 'geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_findmatchingcomps'')', ...
                       'title', 'pop_findmatchingcomps - find matching components marked for rejection', 'mode', 'normal', 'geomvert', [1 1 2 1 1]);
    if isempty(result), return; end
    options = { 'dataset' str2num(result{1}) 'matchcomps' str2num( result{2} ) 'corrthresh' str2double( result{3} ) };
else
    options = varargin;
end

% decode input parameters
% -----------------------
g = finputcheck( options, { 'dataset'       ''         []         [];
                            'corrthresh'    'float'    [0 1]      0.92;
                            'matchcomps'    'float'    []         [];
                            'rejflag'       'float'    [0 1]      [0];
                            'verbose'       'string' { 'on' 'off' } 'on'
                            'nomatcherror'  'string' { 'on' 'off' } 'off' }, 'pop_findmatchingcomps');
if isstr(g), error(g); end

if isempty(g.dataset) && isempty(g.matchcomps)
    error('Either an ALLEG dataset index or set of component scalp maps must be provided as input.');
end

if ~isempty(g.dataset) && ~isempty(g.matchcomps)
    error('Both option can not be provided as input.');
end
outputflag = ~isempty(g.matchcomps); %Checking if input provided for later use
    
% find components in current dataset matching components tagged for rejection 
% in the other dataset
if ~isempty(g.dataset)
    if isstruct(g.dataset)
        TMPEEG = g.dataset;
    else
        global ALLEEG;
        if g.dataset < 1 || g.dataset > length(ALLEEG)
            error('Dataset index out of range');
        end
        TMPEEG = ALLEEG(g.dataset);
    end
    if isempty(TMPEEG.icawinv)
        error(sprintf('Dataset %d does not contain component maps (icawinv)',...
					g.dataset));
    end
    if sum(TMPEEG.reject.gcompreject) == 0
        error( [ 'Matching dataset does not contain components marked for rejection.' 10 ...
               'Follow menu item Tools > Reject data using ICA > Reject components by map' 10 ...
               'to tag non-brain artifact components in that dataset' ] );
    end
        
    g.matchcomps = TMPEEG.icawinv(:,find(TMPEEG.reject.gcompreject));
end

% find matching component
if size(EEG.icawinv,1) ~= size(g.matchcomps,1) && size(EEG.icawinv,1) ~= size(g.matchcomps,2)
    error( [  'The tagged components do not have the same number of' 10 ...
              'channels as those in the current dataset' ]);
end

if size(g.matchcomps,1) ~= size(EEG.icawinv,1) g.matchcomps = g.matchcomps'; end
[corr,indx] = matcorr(EEG.icawinv', g.matchcomps');
corrSelected  = corr(1:size(g.matchcomps,2)) > g.corrthresh;
matchic    = indx(corrSelected);
if outputflag, matchinput = find(corrSelected); else, matchinput = []; end
res = sprintf('Matched %d of %d tagged components using abs. correlation threshold %1.2f', sum(corrSelected), length(corrSelected), g.corrthresh);
if strcmpi(g.verbose, 'on')
    fprintf([ res '\n' ]);
end
if sum(corrSelected) < length(corrSelected) && strcmpi(g.nomatcherror, 'on')
    error([ res 10 'Not all tagged components were matched.' ]);
end
if g.rejflag
    EEG.reject.gcompreject(matchic) = 1;
end

if nargout > 1
    if ~isempty(g.dataset)
        com = sprintf('EEG = pop_findmatchcomps(EEG, ''dataset'', ALLEEG(%d), ''corrthresh'', %1.4f, ''rejflag'',%1.4f);', g.dataset, g.corrthresh,g.rejflag);
    else
        com = sprintf('EEG = pop_findmatchgcomps(EEG, %s);', vararg2str(options));
    end
end
