% std_specplot() - plot STUDY component cluster spectra, either mean spectra 
%                  for all requested clusters in the same figure, with spectra 
%                  for different conditions (if any) plotted in different colors, 
%                  or spectra for each specified cluster in a separate figure 
%                  for each condition,  showing the cluster component spectra plus 
%                  the mean cluster spectrum (in bold). The spectra can be 
%                  plotted only if component spectra have been computed and 
%                  saved with the EEG datasets in Matlab files "[datasetname].icaspec" 
%                  using pop_preclust() or std_preclust(). Called by pop_clustedit(). 
%                  Calls std_readspec() and internal function std_plotcompspec()
% Usage:    
%  >> [STUDY] = std_specplot(STUDY, ALLEEG, key1, val1, key2, val2, ...);  
%  >> [STUDY specdata specfreqs pgroup pcond pinter] = std_specplot(STUDY, ALLEEG, ...);
%
% Inputs:
%   STUDY      - STUDY structure comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - vector of EEG dataset structures for the dataset(s) in the STUDY, 
%                typically created using load_ALLEEG().  
% Optional inputs for component plotting:
%   'clusters' - [numeric vector|'all'] indices of clusters to plot.
%                If no component indices ('comps' below) are given, the average 
%                spectrums of the requested clusters are plotted in the same figure, 
%                with spectrums for different conditions (and groups if any) plotted 
%                in different colors. In 'comps' (below) mode, spectrum for each 
%                specified cluster are plotted in separate figures (one per 
%                condition), each overplotting cluster component spectrum plus the
%                average cluster spectrum in bold. Note this parameter has no effect 
%                if the 'comps' option (below) is used. {default: 'all'}
%   'comps'    - [numeric vector|'all'] indices of the cluster components to plot.
%                Note that 'comps', 'all' is equivalent to 'plotsubjects', 'on'.
%
% Optional inputs for channel plotting:
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel spectrum is plotted (using the 
%                same format as for the cluster component means described above)
%   'subject'  - [numeric vector]  In 'changrp' mode (above), index of 
%                the subject(s) to plot. Else by default, plot all components 
%                in the cluster.
%   'plotsubjects' - ['on'|'off'] When 'on', plot spectrum of all subjects.
%
% Other optional inputs:
%   'plotmode'  - ['normal'|'condensed'] 'normal'  -> plot in a new figure; 
%                 'condensed' -> plot all curves in the current figure in a 
%                 condensed fashion {default: 'normal'}
%   'key','val' - All optional inputs to pop_specparams() are also accepted here
%                 to plot subset of time, statistics etc. The values used by default
%                 are the ones set using pop_specparams() and stored in the
%                 STUDY structure.
% Outputs:
%   STUDY      - the input STUDY set structure with the plotted cluster mean spectra
%                added?? to allow quick replotting.
%   specdata   - [cell] spectral data for each condition, group and subjects.
%                size of cell array is [nconds x ngroups]. Size of each element
%                is [freqs x subjects] for data channels or [freqs x components]
%                for component clusters. This array may be gicen as input 
%                directly to the statcond() function or std_stats() function
%                to compute statistics.
%   specfreqs  - [array] Sprectum point frequency values.
%   pgroup     - [array or cell] p-values group statistics. Output of the 
%                statcond() function.
%   pcond      - [array or cell] condition statistics. Output of the statcond() 
%                function.
%   pinter     - [array or cell] groups x conditions statistics. Output of
%                statcond() function.
%   Example:
%            >> [STUDY] = std_specplot(STUDY,ALLEEG, 'clusters', 2, 'mode', 'apart');
%               % Plot component spectra for STUDY cluster 2, plus the mean cluster 
%               % spectrum (in bold). 
%
%  See also  pop_clustedit(), pop_preclust() std_preclust(), pop_clustedit(), std_readspec()
%
% Authors: Arnaud Delorme, CERCO, August, 2006

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [STUDY, specdata, allfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_specplot;
    return;
end

[STUDY, specdata, allfreqs, pgroup, pcond, pinter] = std_erpplot(STUDY, ALLEEG, 'datatype', 'spec', 'unitx', 'Hz', varargin{:});
