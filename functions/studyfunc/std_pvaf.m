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

function [pvafAve pvafs] = std_pvaf(STUDY, ALLEEG, cluster, varargin);

if nargin < 3
    help std_pvaf;
    return;
end;

[opt addOptions] = finputcheck(varargin, { 'design'  'integer'   []   STUDY.currentdesign; 
                                           'rmclust' 'integer'   []   []; 
                                           'interp'  'struct'    { }   struct([]);
                                           'rmcomps' 'cell'      []    cell(1,length(ALLEEG)) }, 'std_pvaf', 'ignore');
if isstr(opt), error(opt); end;

DES = STUDY.design(opt.design);
for iCell = 1:length(DES.cell)
    [EEG complist rmlistcomp] = std_getdataset(STUDY, ALLEEG, 'cell', iCell, 'cluster', cluster, 'interp', opt.interp, ...
                                    'rmcomps', opt.rmcomps{iCell}, 'rmclust', opt.rmclust, 'interpcomponent', 'on' );
    %EEG = std_getdataset(STUDY, ALLEEG, 'cell', iCell);
    if ~isempty(EEG.icaweights)
         pvafs(iCell) = eeg_pvaf(EEG, [1:size(EEG.icaweights,1)], addOptions{:});
         %pvafs(iCell) = eeg_pvaf(EEG, complist, 'artcomps', rmlistcomp, addOptions{:});
    else pvafs(iCell) = NaN;
    end;
end;
pvafAve = nan_mean(pvafs);

    