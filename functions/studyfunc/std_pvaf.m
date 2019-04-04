% std_pvaf() - Compute 'percent variance accounted for' (pvaf) by specified 
%              ICA component clusters. This function computes eeg_pvaf on each
%              of the component of the cluster and then average them. See 
%              eeg_pvaf for more information. This function uses the
% Usage:
%              >> [pvaf pvafs] = std_pvaf(STUDY, ALLEEG, cluster, 'key', 'val');
% Inputs:
%    EEG       - EEGLAB dataset. Must have icaweights, icasphere, icawinv, icaact.
%    comps     - vector of component indices to sum {default|[] -> progressive mode}
%                In progressive mode, comps is first [1], then [1 2], etc. up to
%                [1:size(EEG.icaweights,2)] (all components); here, the plot shows pvaf.
%
% Optional inputs: 
%   'design'     - [integer] selected design. Default is the current design.
%   'rmcomps'    - [integer array] remove artifactual components (this entry
%                  is ignored when plotting components). This entry contains 
%                  the indices of the components to be removed. Default is none.
%   'interp'     - [struct] channel location structure containing electrode
%                  to interpolate ((this entry is ignored when plotting 
%                  components). Default is no interpolation.
%   Other optional inputs are the same as eeg_pvaf()
%
% Outputs:
%    pvaf      - (real) percent total variance accounted for by the summed 
%                back-projection of the requested clusters.
%    pvafs     - [vector] pvaf for each of the cell of the selected design.
%
% Author:  Arnaud Delorme, SCCN, INC, UCSD, 2012-

% Copyright (C) 2012 Arnaud Delorme, SCCN, INC, UCSD
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

function [pvafAve pvafs] = std_pvaf(STUDY, ALLEEG, cluster, varargin);

if nargin < 3
    help std_pvaf;
    return;
end

[opt addOptions] = finputcheck(varargin, { 'design'  'integer'   []   STUDY.currentdesign; 
                                           'rmclust' 'integer'   []   []; 
                                           'interp'  'struct'    { }   struct([]);
                                           'rmcomps' 'cell'      []    cell(1,length(ALLEEG)) }, 'std_pvaf', 'ignore');
if ischar(opt), error(opt); end

DES = STUDY.design(opt.design);
for iCell = 1:length(DES.cell)
    [EEG complist rmlistcomp] = std_getdataset(STUDY, ALLEEG, 'cell', iCell, 'cluster', cluster, 'interp', opt.interp, ...
                                    'rmcomps', opt.rmcomps{iCell}, 'rmclust', opt.rmclust, 'interpcomponent', 'on' );
    %EEG = std_getdataset(STUDY, ALLEEG, 'cell', iCell);
    if ~isempty(EEG.icaweights)
         pvafs(iCell) = eeg_pvaf(EEG, [1:size(EEG.icaweights,1)], addOptions{:});
         %pvafs(iCell) = eeg_pvaf(EEG, complist, 'artcomps', rmlistcomp, addOptions{:});
    else pvafs(iCell) = NaN;
    end
end
pvafAve = nan_mean(pvafs);

    
