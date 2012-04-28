% STUDY = pop_std_envtopo() - called by eegplugin_std_envtopo. Runs
%                             std_envtopo. Checks if
%                             STUDY.cluster(n).erpdata is present, and if
%                             not all of clusters have it it runs
%                             std_readerp to load and store them.
%
% Usage:
%   >>  STUDY = pop_std_envtopo(STUDY, ALLEEG, 'key', val);
%
% Inputs:
%   STUDY   - an EEGLAB STUDY structure containing EEG structures
%   ALLEEG  - the ALLEEG data structure; can also be an EEG dataset structure.
%
% Optional inputs:
%   same as std_envtopo
%
% Outputs:
%   STUDY  - an EEGLAB STUDY structure containing EEG structures
%
% Author: Makoto Miyakoshi, Arnaud Delorme JSPS/SCCN, INC, UCSD
%
% See also: eegplugin_std_envtopo, pop_envtopo, std_envtopo, envtopo

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2011, Makoto Miyakoshi
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
%
% History
% 04/25/2012 ver 1.9 by Makoto. Description of how to apply low-pass filter added. History function disabled being unable to solve a bug.
% 01/24/2012 ver 1.8 by Makoto. Anti-crash for the cancel bottun. eegh support being attempted
% 01/17/2012 ver 1.7 by Makoto. Input box for 'diff' added.
% 01/13/2012 ver 1.6 by Makoto. Limitation for reduced groups added.
% 01/12/2012 ver 1.5 by Makoto. Limitation is changed to CAUTION.
% 12/26/2011 ver 1.4 by Makoto. Limitation notices added.
% 10/11/2011 ver 1.3 by Makoto. Line75- fixed.
% 10/10/2011 ver 1.2 by Makoto. Two types of scalp topos can be presented.
% 10/06/2011 ver 1.1 by Makoto. Line 77 'design', STUDY.currentdesign.
% 08/01/2011 ver 1.0 by Makoto. Created.

function [STUDY com] = pop_std_envtopo(STUDY, ALLEEG, varargin);

com = '';
if nargin < 2
    help pop_std_envtopo;
    return;
end;

if nargin < 3
    result = inputgui('title', 'Plot envtopo', 'geom', ...
        { {2 12 [0 0] [1 1]} {2 12 [1 0] [1 1]} ...
        {2 12 [0 1] [1 1]} {2 12 [1 1] [1 1]} ...
        {2 12 [0 2] [1 1]} {2 12 [1 2] [1 1]} ...
        {2 12 [0 3] [1 1]} {2 12 [1 3] [1 1]} ...
        {2 12 [0 4] [1 1]} {2 12 [1 4] [1 1]} ...
        {2 12 [0 5] [1 1]} {2 12 [1 5] [1 1]} ...
        {1 12 [0 6] [1 1]}...
        {2 12 [0 7] [1 1]} {2 12 [1 7] [1 1]} ...
        {2 12 [0 8] [1 1]} {2 12 [1 8] [1 1]} ...
        {1 12 [0 10] [1 1]}...
        {1 12 [0 11] [1 1]}...
        {1 12 [0 12] [1 1]}...
        }, 'uilist', ...
        { { 'style' 'text' 'string' 'Time range (in ms) to plot (Ex: -200 500, not -200:500):' } { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' 'Time range (in ms) to rank cluster contributions:' } { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' 'Number of largest contributing clusters to plot (Default: 7)' } { 'style' 'edit' 'string' '7' } ...
        { 'style' 'text' 'string' 'Else plot these cluster numbers only (Ex: 2:4,7):' } { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' 'Cluster numbers to remove from data before plotting:' } { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' 'Difference. 1 1 2 1 = Cond1 Group1 minus Cond2 Group1.'} { 'style' 'edit' 'string' '' }...
        { 'style' 'text' 'string' '        (Provide 4 numbers. If no Condition/Group, enter 1.)'}...
        { 'style' 'text' 'string' 'Include only the clusters that were part of the clustering' } { 'style' 'checkbox', 'value', 1} ...
        { 'style' 'text' 'string' 'Optional inputs (''amplimits'', ''fillclust'', etc. see std_envtopo.m)' } { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' '++++++++++++++++  The range options in ''Params'' for ERP MUST be empty. Low-pass filter is available from'}...
        { 'style' 'text' 'string' '++++ CAUTION!!! ++++   ''Params'' for ERP. Group MUST be a second variable in design. DO NOT exclude any group' }...
        { 'style' 'text' 'string' '++++++++++++++++   in the design interface. DO NOT change the name ''outlier 2''. Clusters MUST NOT be empty.'}});
    
    % if canceled, escape
    if isempty(result), return, end
    
    % prepare optional inputs for std_envtopo
    options = '';
    if ~isempty( result{1} ), options = [ options '''timerange'',[' result{1} '],' ]; end;
    if ~isempty( result{2} ), options = [ options '''limcontrib'',[' result{2} '],' ]; end;
    if ~isempty( result{3} ), options = [ options '''clustnums'',[-' result{3} '],' ]; end;
    if ~isempty( result{4} ), options = [ options '''clustnums'',[' result{4} '],' ]; end;
    if ~isempty( result{5} ), options = [ options '''clust_exclude'',[' result{5} '],' ]; end;
    if ~isempty( result{6} ), options = [ options '''diff'',[' result{6} '],' ]; end;
    if result{7} == 1; options = [ options '''onlyclus'', ''on''']; else options = [ options '''onlyclus'', ''off''']; end
    if ~isempty( result{8} ), options = [ options ',' result{8} ]; end;
    arguments = eval([ '{' options '}' ]);
else
    arguments = varargin;
end;

if length(STUDY.cluster) == 1
    errordlg2('Cannot plot envtopo with the parent cluster only');
end;

for n = 2:length(STUDY.cluster)
    [STUDY, datavals, xvals, setinds, allinds] = std_readerp(STUDY, ALLEEG, 'design', STUDY.currentdesign, 'clusters', n, 'singletrials', 'off');
end

% run std_evntopo
std_envtopo(STUDY, ALLEEG, arguments{:});
% com = ['std_envtopo(STUDY, ALLEEG,' options ');'];
