% std_envtopo() - Creates an envtopo() image for a STUDY set using component cluster
%                 contributions instead of individual components.  Plots the envelope
%                 of the data epoch grand mean ERP, plus envelopes and average scalp maps
%                 for specified or largest-contributing clusters for each condition.
%                 Click on individual axes to examine them in detail (using axcopy()).
%                 See envtopo() for further details.
%
% Usage:          
%                >> std_envtopo(STUDY, ALLEEG, 'key1', 'val1', ...);
%
% Inputs:
%   STUDY        = an EEGLAB STUDY structure containing EEG structures
%   ALLEEG       = the ALLEEG data structure; can also be an EEG dataset structure.
%
% Optional inputs:
%  'clustnums'   = [integer array] vector of cluster numbers to plot.  Else if
%                   int < 0, the number of largest contributing clusters to plot
%                   {default|[] -> 7}
%  'conditions'  = [integer array] vector of condition indices to plot
%  'subclus'     = [integer array] vector of cluster numbers to omit when  computing
%                    the ERP envelope of the data (e.g., artifact
%                    clusters). By default, these clusters are excluded from
%                   the grandERP envelope. See the next option 'onlyclus'.
%  'onlyclus'    = [ 'on' | 'off'] dataset components to include in the grand ERP.
%                  'on' will include only the components that were part of the
%                   clustering. For example, if components were rejected from
%                   clustering because of high dipole model residual variance,
%                   don't include their data in the grand ERP.
%                  'off' will include all components in the datasets except
%                   those in the subtructed ('subclus') clusters {default 'on'}.
%  'baseline'    = [minms maxms] - a new baseline to remove from the grand
%                   and cluster ERPs.
%  'diff'        = [condition1 group 1 condition 2 group 2] Perform subtractiono
%                   between conbination of condition/group. Ex. 'diff', [1 1 2 1]
%                   performs condition1_group1 minus condition2_group1.
%                   Must be provided with 4 numbers. If no condition/group, 1.
%  'timerange'   = data epoch start and end input latencies (in ms)
%                   {default: from 'limits' if any}
%  'amplimits'      = [minuV maxuV]. {default: use data uV limits}
%  'limcontrib'  = [minms maxms]  time range (in ms) in which to rank cluster contributions
%                   (boundaries = thin dotted lines) {default|[]|[0 0] -> plotting limits}
%  'vert'        = vector of times (in ms) at which to plot vertical dashed lines
%                   {default|[] -> none}
%  'sortvar'     = ['pvaf'|'rv'] if 'pvaf', sort by percent variance accounted for.
%                   If 'rv', sort by relative variance.
%                   pvaf(component) = 100-100*variance(data-cluster))/variance(data)
%                   rv(component)   = 100*variance(component)/variance(data) {default: 'rv'}
%  'topotype'    = ['inst'|'ave'] If 'inst', show the instantaneous map at the specific timepoint specified by the line.
%                   If 'ave', show that of the mean which is the same as the
%                   stored topomap. {default: 'inst'}
%  'fillclust'    = [integer] fill the numbered component envelope with red. {default|[]|0 -> no fill}
%
% Author: Makoto Miyakoshi, Hilit Serby, Arnold Delorme, Scott Makeig
% The original version of this script was written by Hilit Serby and redesigned by Makoto Miyakoshi.
%
% See also: eegplugin_std_envtopo, pop_std_envtopo, std_envtopo, envtopo

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2011, Makoto Miyakoshi, Hilit Serby, Arnold Delorme, Scott Makeig
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA
%
% History
% 04/29/2012 ver 5.2 by Makoto. A bug fixed. STUDY.design(STUDY.currentdesign)
% 04/24/2012 ver 5.1 by Makoto. Revived 'amplimit' 'fillclust' options, and vertical doted lines for limcontrib range. Default changed into 'pvaf' from 'rv'. 
% 02/28/2012 ver 5.0 by Makoto. Values are now normalized by the number of subjects participating to the cluster. Group difference is better addressed when performing subtraction across groups. A group ratio is output too.
% 02/24/2012 ver 4.5 by Makoto. New color schema to improve color readability & shared line colors across multiple plots
% 02/22/2012 ver 4.4 by Makoto. Now plot title can accept numbers
% 01/24/2012 ver 4.3 by Makoto. Output added for the function.
% 01/17/2012 ver 4.2 by Makoto. Improved diff option; now mutually exclusive with normal plot. Diff title improved too. 
% 01/06/2012 ver 4.1 by Makoto. Time range options fixed. STUDY.design selection fixed. Wrong labels fixed. (Thanks to Agatha Lenartowicz!) 'sortvar' and 'topotype' implemented in the main function. Help edited.
% 12/30/2011 ver 4.0 by Makoto. Main function 100% redesigned. Full support for Study design.
% 12/26/2011 ver 3.2 by Makoto. Formats corrected.
% 12/25/2011 ver 3.1 by Makoto. Bug fixed about choosing wrong clusters under some conditions. 
% 10/12/2011 ver 3.0 by Makoto. Completed the alpha version. 
% 10/11/2011 ver 2.2 by Makoto. Time range, group with no entry fixed
% 10/10/2011 ver 2.1 by Makoto. Two types of scalp topos can be presented.
% 10/08/2011 ver 2.0 by Makoto. Multiple groups supported (but one group one time using option). Result topos are now retrieved from the stored ones in STUDY 
% 08/01/2011 ver 1.0 by Makoto. Updated the original script to read data from memory.

function STUDY = std_envtopo(STUDY, ALLEEG, varargin);

% if there is < 2 arguments, show help
if nargin < 2
    help std_envtopo;
    return
end

% if arguments are in odd number, wrong input
if mod(nargin,2) % if not an even number of arguments
    error('std_envtopo: Input argument list must be pairs of: ''keyx'', ''valx'' ');
end

% run finputcheck
arglist = finputcheck(varargin, ...
             {'clustnums'     'integer'  [-length(STUDY.cluster):length(STUDY.cluster)  ]     [];...
              'conditions'    'integer'  [1:length(STUDY.design(STUDY.currentdesign).variable(1,1).value)] [1:length(STUDY.design(STUDY.currentdesign).variable(1,1).value)];...
              'clust_exclude' 'integer'  [1:length(STUDY.cluster)]                            [];...
              'onlyclus'      'string'   {'on', 'off'}                                        [];...
              'baseline'      'real'     []                                                   [];...
              'diff'          'integer'  []                                                   [];...
              'timerange'     'real'     []               [ALLEEG(1,1).times([1 end])];...
              'amplimits'     'real'     []                                                   [];...
              'limcontrib'    'real'     []                                                   [];...
              'vert'          'real'     []                                                   [];...
              'sortvar'       'string'   {'pvaf', 'rv'}                                   'pvaf';...
              'topotype'      'string'   {'inst', 'ave'}                                  'inst';...
              'fillclust'      'integer'  []                                                   []});
          
if isstr(arglist)
    error(arglist);
end
clear varargin

% Timepoints to plot (from sample)
arglist.timepoints = find(STUDY.cluster(1,2).erptimes >= arglist.timerange(1) & STUDY.cluster(1,2).erptimes <= arglist.timerange(2));

% Convert limcontrib period to timepoints (from sample)
if ~isempty(arglist.limcontrib)
    arglist.envtimelimits = find(STUDY.cluster(1,2).erptimes >= arglist.limcontrib(1) & STUDY.cluster(1,2).erptimes <= arglist.limcontrib(2));
end

% Convert baseline period to timepoints (from sample)
if  ~isempty(arglist.baseline)
    arglist.baselinelimits = find(STUDY.cluster(1,2).erptimes >= arglist.baseline(1) & STUDY.cluster(1,2).erptimes <= arglist.baseline(2));
end

% Outlier cluster & all clusters included
if strcmp(STUDY.cluster(2).name(1:5), 'outli')
    arglist.clust_all = [3:length(STUDY.cluster)];% exclude parent (1) and outliers (2)
else
    arglist.clust_all = [2:length(STUDY.cluster)];% exclude parent (1) only
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Determine clusters to plot    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if arglist.clustnums > 0;
    arglist.clust_selected = arglist.clustnums;
    arglist.clust_topentries = [];
elseif arglist.clustnums < 0
    arglist.clust_selected = [];
    arglist.clust_topentries = abs(arglist.clustnums);
end
if strcmp(arglist.onlyclus, 'on')
    arglist.clust_grandERP = setdiff(arglist.clust_all, arglist.clust_exclude);
else
    arglist.clust_grandERP = arglist.clust_all;
end
arglist.clustlabels = {STUDY.cluster(1,arglist.clust_grandERP).name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Separate topos into groups     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(STUDY.group) > 1
    [STUDY, ALLEEG] = group_topo_now(STUDY, ALLEEG);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Convolve scalp topo with ERP   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

for cls = 1:length(arglist.clust_grandERP) % For all clusters for grandERP
    for columnwise = 1:size(STUDY.cluster(1,arglist.clust_grandERP(cls)).erpdata, 1) % For the first variable
        for rowwise = 1:size(STUDY.cluster(1,arglist.clust_grandERP(cls)).erpdata, 2) % For the second variable
            if ~isempty(STUDY.cluster(1,arglist.clust_grandERP(cls)).erpdata{columnwise, rowwise}) % Detect no group entry
                erp   = STUDY.cluster(1,arglist.clust_grandERP(cls)).erpdata{columnwise, rowwise}';
                if  length(STUDY.group) > 1
                    topo  = STUDY.cluster(1,arglist.clust_grandERP(cls)).topoall_group{1, rowwise};
                else
                    topo  = STUDY.cluster(1,arglist.clust_grandERP(cls)).topoall;
                end
                for n = 1:length(topo)
                    tmp1 = reshape(topo{1,n}, size(topo{1,n},1)*size(topo{1,n},2), 1);
                    tmp1 = tmp1(find(~isnan(tmp1)));
                    tmp2(:,n) = tmp1;
                end
                topo = tmp2; clear n tmp1 tmp2
                topoerpconv{1, cls}{columnwise,rowwise} = topo*erp; % Convolve topo with erp
            else
                disp([char(10) ' CAUTION: No IC entry in Cluster ' num2str(arglist.clust_grandERP(1,1)+cls-1) ' Group ' num2str(rowwise)])
                topoerpconv{1, cls}{columnwise,rowwise} = zeros(3409, length(arglist.timepoints)); % Dummy data
            end
            topoerpconv_pile(:,:,cls, columnwise, rowwise) = topoerpconv{1, cls}{columnwise,rowwise};
        end
    end
end
clear cls columnwise erp rowwise topo topoerpconv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Normalization by dividing with the number of subjects %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setinds = {STUDY.cluster.setinds};
for n = 1:length(setinds)
    for row = 1:size(setinds{1,1}, 1)
        for column = 1:size(setinds{1,1}, 2)
            nsubject{1,n}(row,column) = length(unique(setinds{1,n}{row,column}));
        end
    end
    nsubject_group{1,n} = nsubject{1,n}(1,:);
end
arglist.normalize_cluster = nsubject_group;

for n = 1:size(topoerpconv_pile, 3) % for the number of clusters
        str = sprintf('In %s, ', arglist.clustlabels{1,n});
    for m = 1:size(setinds{1,1}, 2) % for the number of groups
        topoerpconv_pile(:,:,n,:,m) = topoerpconv_pile(:,:,n,:,m)./arglist.normalize_cluster{1,arglist.clust_grandERP(n)}(1,m);
        str = [str sprintf('%d %s ', arglist.normalize_cluster{1,arglist.clust_grandERP(n)}(1,m), num2str(STUDY.group{1,m}))];
    end
    if length(STUDY.group) > 1
        disp(str)
    end
end
clear column n row setinds nsubject* str

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply new baseline (if specified) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

if isfield(arglist, 'baselinelimits')
    topoerpconv_pile = topoerpconv_pile - repmat(mean(topoerpconv_pile(:, arglist.baselinelimits,:,:,:), 2), [1, size(topoerpconv_pile, 2), 1,1,1]);
    fprintf('\nNew baseline applied.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate outermost envelope for each conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

outermostenv = sum(topoerpconv_pile, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select clusters (if specified) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if arglist.clustnums > 0;
    [dummy indx] = intersect(arglist.clust_grandERP, arglist.clustnums);
    topoerpconv_pile = topoerpconv_pile(:,:,indx,:,:);
    tmp = {arglist.clustlabels(1, indx)};
    arglist.clustlabels = tmp{1,1}; clear tmp dummy indx
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Plot envtopo for each condition & group    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% Create an dummy color log to share line colors for the same clusters across plots
colorlog = zeros(1,2);

if isempty(arglist.diff)
    for columnwise = 1:size(topoerpconv_pile, 4) % For all conditions
        for rowwise = 1:size(topoerpconv_pile, 5) % For all groups
            figure; set(gcf,'Color', [0.93 0.96 1]); orient landscape;
            arglist.title = [num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,1).value{1, columnwise}) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).value{1, rowwise})];
            colorlog = envtopo_plot(STUDY, ALLEEG, squeeze(outermostenv(:,arglist.timepoints,1,columnwise,rowwise)), squeeze(topoerpconv_pile(:,arglist.timepoints,:,columnwise,rowwise)), colorlog, arglist);
        end
    end
end
clear columnwise rowwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Else, plot difference envtopo    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

if ~isempty(arglist.diff)
    % First, plot title is determined
    if     length(STUDY.design(STUDY.currentdesign).variable(1,1).value) == 1 && arglist.diff(1) == 1 && arglist.diff(3) == 1 % Group only if no condition
        arglist.title = [num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.diff(1,2)}) ' minus ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.diff(1,4)})];
    elseif     length(STUDY.group) == 1 && arglist.diff(2) == 1 && arglist.diff(4) == 1 % Condition only if no group
        arglist.title = [num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.diff(1,1)}) ' minus ' num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.diff(1,3)})];
    else   % if both condition and group
        arglist.title = [num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.diff(1,1)}) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.diff(1,2)}) ' minus ' num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.diff(1,3)}) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.diff(1,4)})];
    end
    % Calculation and plotting
    figure; set(gcf,'Color', [0.93 0.96 1]); orient landscape;
    outermostenv = outermostenv(:,:,:,arglist.diff(1,1),arglist.diff(1,2)) - outermostenv(:,:,:,arglist.diff(1,3),arglist.diff(1,4));
    topoerpconv_pile = topoerpconv_pile(:,:,:,arglist.diff(1,1),arglist.diff(1,2)) - topoerpconv_pile(:,:,:,arglist.diff(1,3),arglist.diff(1,4));
    envtopo_plot(STUDY, ALLEEG, squeeze(outermostenv(:,arglist.timepoints,:,:,:)), squeeze(topoerpconv_pile(:,arglist.timepoints,:,:,:)), colorlog, arglist);
end

% Pass the current STUDY to the base workspece (this is tricky...)
userdat = get(2, 'UserData');
userdat{1}{2} = STUDY;
set(2, 'UserData', userdat);

% envtopo_plot() - Plot the envelope of a data epoch, plus envelopes and scalp maps of specified
%             or largest-contributing components. If a 3-D input matrix, operates on the
%             mean of the data epochs. Click on individual axes to examine them in detail.
% Usage:
%             >> envtopo_plot(grandenvERP,projERP);
%             >> [compvarorder,compvars,compframes,comptimes,compsplotted,pvaf] ...
%                                           = envtopo(grandERP,projERP, 'key1', val1, ...);
% Inputs:
%  grandERP     = The grand average ERP (chans, frames),
%                              (see std_granderp() and std_envtopo() for details).
%  projERP      = A cell array of the projected ERPs of the desired components, each cell size is (chans,frames),
%                              (see std_granderp() and std_envtopo() for details).
%
% Optional inputs:
%  'clustnums'  = [integer array] vector of clusters numbers to plot {default|0 -> all}
%                  Else if int < 0, the number of largest contributing clusters to plot
%                  {default|[] -> 7}
%  'timerange' = start and end input data latencies (in ms) {default: from 'limits' if any}
%  'amplimits'    = [minuV maxuV].
%  'limcontrib' = [minms maxms]  time range (in ms) in which to rank component contribution
%                  (boundaries shown with thin dotted lines) {default|[]|[0 0] -> plotting limits}
%  'title'      = [string] plot title {default|[] -> none}
%  'plotchans'  = [integer array] data channels to use in computing contributions and envelopes,
%                  and also for making scalp topo plots {default|[] -> all}
%  'voffsets'   = [float array] vertical line extentions above the data max to disentangle
%                  plot lines (left->right heads, values in y-axis units) {def|[] -> none}
%  'colors'     = [string] filename of file containing colors for envelopes, 3 chars
%                  per line, (. = blank). First color should be "w.." (white)
%                  Else, 'bold' -> plot default colors in thick lines.
%                  {default|[] -> standard Matlab color order}
%  'vert'       = vector of times (in ms) at which to plot vertical dashed lines {default|[] -> none}
%  'icawinv'    = [float array] inverse weight matrix. By default computed by inverting
%                  the weight matrix (but if some components have been removed, then
%                  weight's pseudo-inverse matrix does not represent component's maps).
%  'icaact'     = [float array] component activations. By default these are computed
%                  from the input weight matrix.
%  'envmode'    = ['avg'|'rms'] compute the average envelope or the root mean square
%                  envelope {default: 'avg'}
%  'subcomps'   = [integer vector] indices of components to remove from data before
%                  plotting. {default: none}
%  'clustlabels'   = [cell array of strings] the size of the clusters number, to label
%                  the clusters with labels. {default: none}
%  'sumenv'     = ['on'|'off'|'fill'] 'fill' -> show the filled envelope of the summed projections
%                  of the selected components; 'on' -> show the envelope only {default: 'fill'}
%  'actscale'   = ['on'|'off'] scale component scalp maps by maximum component activity in the
%                  designated (limcontrib) interval. 'off' -> scale scalp maps individually using
%                  +/-max(abs(map value)) {default: 'off'}
%  'dispmaps'   = ['on'|'off'] display component numbers and scalp maps {default: 'on'}
%  'topoplotkey','val' = any optional arguments for topoplot.
%
% Outputs:
%  compvarorder = component numbers in decreasing order of max variance in data
%  compvars     = component max variances
%  compframes   = frames of max variance
%  comptimes    = times of max variance
%  compsplotted = components plotted
%  mv|pvaf|rv   = max variance, percent variance accounted for, or relative variance (see 'sortvar')
%
% Notes:
%  To label maps with other than component numbers, put four-char strings into a local (pwd) file
%  named 'envtopo.labels' (using . = space) in time-order of their projection maxima
%
% Authors: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 3/1998
%          Edited for std_envtopo by Makoto Miyakoshi.
%
% See also: timtopo()

% Copyright (C) 3-10-98 from timtopo.m Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

function [colorlog] = envtopo_plot(STUDY, ALLEEG, grandERP,projERP, colorlog, varargin)

% data format conversion for compatibility
g = varargin{1,1}; clear varargin
for n = 1:size(projERP, 3);
    tmp{1,n} = projERP(:,:,n);
end
projERP = tmp;
clear n tmp

% for backwards compatibility 11/2004 -sm
all_bold = 0;
BOLD_COLORS = 1;    % 1 = use solid lines for first 5 components plotted
% 0 = use std lines according to component rank only
FILL_COMP_ENV = 0;  % default no fill
FILLCOLOR   = [.815 .94 1]; % use lighter blue for better env visibility
MAXTOPOS = 20;      % max topoplots to plot
VERTWEIGHT = 2.0;  % lineweight of specified vertical lines
LIMCONTRIBWEIGHT = 1.2; % lineweight of limonctrib vertical lines

myfig =gcf;         % remember the current figure (for Matlab 7.0.0 bug)
xmin = 0; xmax = 0;

% add options (default values)
g.envmode     = 'avg';
g.dispmaps    = 'on';
g.voffsets    = [];
g.colorfile   = '';
g.colors      = '';
g.xlabel      = 'on';
g.ylabel      = 'on';
g.actscale    = 'off';
g.topoarg     = 0;
g.sumenv      = 'fill';

% empty cluster crashes the code, so 
for n = 2:length(STUDY.cluster)
    if ~isempty(STUDY.cluster(1,n).topoall)
        g.gridind = find(~isnan(STUDY.cluster(1,n).topoall{1,1}));
        break
    end
end
clear n
%
% Check input flags and arguments
%
[tmp,frames] = size(grandERP);

if ~isempty(g.colors)
    g.colorfile = g.colors; % retain old usage 'colorfile' for 'colors' -sm 4/04
end

if ~isempty(g.vert)
    g.vert = g.vert/1000; % convert from ms to s
end
%
%%%%%% Collect information about the gca, then delete it %%%%%%%%%%%%%
%
uraxes = gca; % the original figure or subplot axes
pos=get(uraxes,'Position');
axcolor = get(uraxes,'Color');
delete(gca)

%
%%% Convert g.timerange and g.limcontrib to sec from ms %%%%
%
g.timerange = g.timerange/1000;   % the time range of the input data
g.limcontrib = g.limcontrib/1000; % the time range in which to select largest components


%
%%%%%%%%%%%% Collect time range information %%%%%%%%%%%%%%%%%%%%%%%%%%
%

xunitframes = 0;
xmin = g.timerange(1);
xmax = g.timerange(2);
pmin = xmin;
pmax = xmax;

dt = (xmax-xmin)/(frames-1);  % sampling interval in sec
times=xmin*ones(1,frames)+dt*(0:frames-1); % time points in sec

%
%%%%%%%%%%%%%%% Find limits of the component selection window %%%%%%%%%
%
if any(g.limcontrib ~= 0)
    if xunitframes
        g.limcontrib = g.limcontrib*1000; % if no time limits, interpret
    end                                  % limcontrib as frames
    if g.limcontrib(1)<xmin
        g.limcontrib(1) = xmin;
    end
    if g.limcontrib(2)>xmax
        g.limcontrib(2) = xmax;
    end
    srate = (frames-1)/(xmax-xmin);
    limframe1  = round((g.limcontrib(1)-xmin)*srate)+1;
    limframe2  = round((g.limcontrib(2)-xmin)*srate)+1;
    g.vert(end+1) =  g.limcontrib(1);
    g.vert(end+1) =  g.limcontrib(2);
else
    limframe1 = 1;
    limframe2 = frames;
end;

%
%%%%%%%%%%%%%%%%%%%%% Line color information %%%%%%%%%%%%%%%%%%%%%
%
% 16 colors names officially supported by W3C specification for HTML
colors{1,1}  = [1 1 1];            % White
colors{2,1}  = [1 1 0];            % Yellow
colors{3,1}  = [1 0 1];            % Fuchsia
colors{4,1}  = [1 0 0];            % Red
colors{5,1}  = [0.75  0.75  0.75]; % Silver
colors{6,1}  = [0.5 0.5 0.5];      % Gray
colors{7,1}  = [0.5 0.5 0];        % Olive
colors{8,1}  = [0.5 0 0.5];        % Purple
colors{9,1}  = [0.5 0 0];          % Maroon
colors{10,1} = [0 1 1];            % Aqua
colors{11,1} = [0 1 0];            % Lime
colors{12,1} = [0 0.5 0.5];        % Teal
colors{13,1} = [0 0.5 0];          % Green
colors{14,1} = [0 0 1];            % Blue
colors{15,1} = [0 0 0.5];          % Navy
colors{16,1} = [0 0 0];            % Black

% Silver is twice brighter because used for background
colors{5,1} = [0.875 0.875 0.875];

% Choosing and sorting 12 colors for line plot, namely Red, Blue, Green, Fuchsia, Lime, Aqua, Maroon, Olive, Purple, Teal, Navy, and Gray
linecolors = colors([4 13 14 3 11 10 9 7 8 12 15 6]);

%
%%%%%%%%%%%%%%%% Check other input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

[tmp,pframes] = size(projERP{1});
if frames ~= pframes
    error('Size of trial in projERP and grandenvERP do not agree');
end

if isempty(g.voffsets) | ( size(g.voffsets) == [1,1] & g.voffsets(1) == 0 )
    g.voffsets = zeros(1,MAXTOPOS);
end

if isempty(g.clustnums) | g.clustnums(1) == 0
   g.clustnums = 1:length(projERP); % by default the number of projected ERP input
end
if min(g.clustnums) < 0
    if length(g.clustnums) > 1
        error('Negative clustnums must be a single integer.');
    end
    if -g.clustnums > MAXTOPOS
        fprintf('Can only plot a maximum of %d components.\n',MAXTOPOS);
        return
    else
        MAXTOPOS = -g.clustnums;
        g.clustnums = 1:length(projERP); % by default the number of projected ERP input
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Process components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ncomps = length(g.clustnums);

%
%%%%%%%%%%%%%%% Compute plotframes and envdata %%%%%%%%%%%%%%%%%%%%%
%
ntopos = length(g.clustnums);
if ntopos > MAXTOPOS
    ntopos = MAXTOPOS; % limit the number of topoplots to display
end

plotframes = ones(ncomps,1);
%
% first, plot the data envelope
%
envdata = zeros(2,frames*(ncomps+1));
envdata(:,1:frames) = envelope(grandERP, g.envmode);

fprintf('Data epoch is from %.0f ms to %.0f ms.\n',1000*xmin,1000*xmax);
fprintf('Plotting data from %.0f ms to %.0f ms.\n',1000*xmin,1000*xmax);
fprintf('Comparing maximum projections for components:  ');
if ncomps>32
    fprintf('\n');
end
compvars = zeros(1,ncomps);
mapsigns = zeros(1,ncomps);

%
% Compute frames to plot
%
sampint  = (xmax-xmin)/(frames-1);     % sampling interval in sec
times    = xmin:sampint:xmax;          % make vector of data time values

[v minf] = min(abs(times-pmin));
[v maxf] = min(abs(times-pmax));
pframes  = minf:maxf;         % frames to plot
ptimes   = times(pframes);    % times to plot
if limframe1 < minf
    limframe1 = minf;
end
if limframe2 > maxf
    limframe2 = maxf;
end

%
%%%%%%%%%%%%%% find max variances and their frame indices %%%%%%%%%%%
%
nvals = size(grandERP,1)*length(limframe1:limframe2);
for c = 1:ncomps
    if ~rem(c,5)
        fprintf('%d ... ',g.clustnums(c)); % c is index into clustnums
    end
    if ~rem(c,100)
        fprintf('\n');
    end

    envdata(:,c*frames+1:(c+1)*frames) = envelope(projERP{c}, g.envmode);

    [maxval,maxi] = max(sum(projERP{c}(:,limframe1:limframe2).*projERP{c}(:,limframe1:limframe2)));
    % find point of max variance for comp c
    compvars(c)   = maxval;
    maxi = maxi+limframe1-1;
    plotframes(c) = maxi;
    maxproj(:,c)  = projERP{c}(:,maxi); % Note: maxproj contains only g.plotchans -sm 11/04

end % component c
fprintf('\n');

%
%%%%%%%%%%%%%%% Compute component selection criterion %%%%%%%%%%%%%%%%%%%%%%%%%%
%

% compute pvaf
if ~xunitframes
    fprintf('  in the interval %3.0f ms to %3.0f ms.\n',1000*times(limframe1),1000*times(limframe2));
end

vardat = var(reshape(grandERP(:,limframe1:limframe2),1,nvals)); % find full data variance in interval
for c = 1:length(projERP)
    % now calculate the pvaf over the whole time period for printing
    if strcmpi(g.sortvar, 'pvaf') | strcmpi(g.sortvar,'pv')
        diffdat = grandERP(:,limframe1:limframe2)-projERP{c}(:,limframe1:limframe2);
        diffdat = reshape(diffdat,1,nvals);
        pvaf(c) = 100-100*(var(diffdat)/vardat); %var of diff div by var of full data
        ot   = 'pvaf';
    elseif strcmpi(g.sortvar, 'rv')% var(clust)/var(data)
        pvaf(c) = 100*(var(reshape(projERP{c}(:,limframe1:limframe2),1,nvals))/vardat);
        ot   = 'rv';
    end;
end

%
%%%%%%%%%%%%%%%%%%%%%%%%% Sort by max variance in data %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[pvaf,compx] = sort(pvaf);  % sort clustnums on max pvaf
pvaf = pvaf(ncomps:-1:1);  % reverse order of sort (max:min)       
compx        = compx(ncomps:-1:1);    % reverse order of sort
compvarorder = g.clustnums(compx);    % actual cluster numbers (output var)
plotframes   = plotframes(compx);     % plotted comps have these max frames

maxproj    = maxproj(:,compx); % maps in plotting order 
compvars = compvars(compx);
if ~isempty(g.clustlabels) 
   complabels = g.clustlabels(compx);  
   complabels = complabels(1:ntopos);% actual component numbers (output var)
end;

%
%%%%%%%%%%%%%%%%%%%%%%%% Reduce to ntopos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
compx = compx(1:ntopos);
pvaf = pvaf(1:ntopos);
compvars = compvars(1:ntopos);% max variances
maxproj    = maxproj(:,1:ntopos);
compsplotted = compvarorder(1:ntopos);% (output var)
compvarorder = compvarorder(1:ntopos);
[plotframes,ifx] = sort(plotframes(1:ntopos));% sort plotframes on their temporal order (min:max)
plottimes  = times(plotframes);       % convert to times in ms
compx      = compx(ifx);              % indices into clustnums, in plotting order
maporder   = compvarorder(ifx);       % reorder cluster numbers
maxproj    = maxproj(:,ifx);          % maps in plotting order
if ~isempty(g.clustlabels) 
   complabels = complabels(ifx);     % actual component numbers (output var)
end;
pvaf = pvaf(ifx);
vlen = length(g.voffsets); % extend voffsets if necessary
while vlen< ntopos
    g.voffsets = [g.voffsets g.voffsets(vlen)]; % repeat last offset given
    vlen=vlen+1;
end

head_sep = 1.2;
topowidth = pos(3)/(ntopos+(ntopos-1)/5); % width of each topoplot
if topowidth > 0.20    % adjust for maximum height
    topowidth = 0.2;
end

if rem(ntopos,2) == 1  % odd number of topos
    topoleft = pos(3)/2 - (floor(ntopos/2)*head_sep + 0.5)*topowidth;
else % even number of topos
    topoleft = pos(3)/2 - (floor(ntopos/2)*head_sep)*topowidth;
end

%
%%%%%%%%%%%%%%%%%%%% Print times and frames of comp maxes %%%%%%%%%%%%%%
%
if ~xunitframes
    fprintf('    with max var at times (ms): ');
    for t=1:ntopos
        fprintf('%4.0f  ',1000*plottimes(t));
    end
    fprintf('\n');
end

fprintf('                  epoch frames: ');
for t=1:ntopos
    fprintf('%4d  ',limframe1-1+plotframes(t));
end
fprintf('\n');

if strcmpi(g.sortvar,'on') | strcmpi(g.sortvar,'pvaf')
    fprintf('    Component pvaf in interval:  ');
    for t=1:ntopos
        fprintf('%s ',num2str(pvaf(t)));
    end
    fprintf('\n');
end

sumproj = zeros(size(projERP{1}));
for n = 1:ntopos
    sumproj = sumproj + projERP{compx(n)}; % add up all cluster projections
end

totlimdat = grandERP(:,limframe1:limframe2);
sumlimdat = sumproj(:,limframe1:limframe2);
diffdat = totlimdat - sumlimdat; 
diffdat = reshape(diffdat,1,nvals);

if strcmpi(g.sortvar,'pvaf')   | strcmpi(g.sortvar,'pv')
    sumpvaf = 100-100*(var(diffdat)/vardat); 
    %sumpvaf = 100-100*(var(reshape(totlimdat-sumlimdat,1,nvals))/vardat); 
    ot   = 'pvaf';
elseif strcmpi(g.sortvar, 'rv')
    sumpvaf = 100*(var(reshape(sumlimdat,1,nvals))/vardat); 
    ot   = 'rv';
end;

if ~xunitframes
    fprintf('    Summed component %s in interval [%s %s] ms: %4.2f%%\n',ot, int2str(1000*times(limframe1)),int2str(1000*times(limframe2)), sumpvaf);
end
%
%%%%%%%%%%%%%%%%%%%%% Plot the data envelopes %%%%%%%%%%%%%%%%%%%%%%%%%
%
BACKCOLOR = [0.7 0.7 0.7];
newaxes=axes('position',pos);
axis off
set(newaxes,'FontSize',16,'FontWeight','Bold','Visible','off');
set(newaxes,'Color',BACKCOLOR); % set the background color
delete(newaxes) %XXX

% site the plot at bottom of the current axes
axe = axes('Position',[pos(1) pos(2) pos(3) 0.6*pos(4)],...
    'FontSize',16,'FontWeight','Bold');

set(axe,'GridLineStyle',':')
set(axe,'Xgrid','off')
set(axe,'Ygrid','on')
axes(axe)
set(axe,'Color',axcolor);

%
%%%%%%%%%%%%%%%%% Plot the envelope of the summed selected components %%%%%%%%%%%%%%%%%
%

if isempty(g.amplimits)
    envx = [1,compx+1];
    ymin = min(matsel(envdata,frames,0,1,envx(1)));
    ymax = max(matsel(envdata,frames,0,2,envx(1)));
else
    ymin = g.amplimits(1);
    ymax = g.amplimits(2);
end

if strcmpi(g.sumenv,'on')  | strcmpi(g.sumenv,'fill')
    sumenv = envelope(sumproj, g.envmode);
    if strcmpi(g.sumenv,'fill')
        %
        % Plot the summed projection filled
        %
        mins = matsel(sumenv,frames,0,2,0);
        p=fill([times times(frames:-1:1)],...
            [matsel(sumenv,frames,0,1,0) mins(frames:-1:1)], colors{5,1});
        set(p,'EdgeColor',colors{5,1});
        hold on
        %
        % Overplot the data envelope so it is not covered by the fill()'d component
        %
        p=plot(times,matsel(envdata,frames,0,1,1),'Color', colors{16,1});% plot the max
        set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)
        p=plot(times,matsel(envdata,frames,0,2,1),'Color', colors{16,1});% plot the min
        set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)

    else % if no 'fill'
        tmp = matsel(sumenv,frames,0,2,0);
        p=plot(times,tmp);% plot the min
        hold on
        set(p,'Color', colors{5,1});
        set(p,'linewidth',2);
        p=plot(times,matsel(sumenv,frames,0,1,0));% plot the max
        set(p,'linewidth',2);
        set(p,'Color', colors{5,1});
    end
end

%
% %%%%%%%%%%%%%%%%%%%%%%%% Plot the computed component envelopes %%%%%%%%%%%%%%%%%%
%
envx = [1,compx+1];
hold on
for c = 1:ntopos+1
    for envside = [1 2]
        if envside == 1 % maximum envelope
            curenv = matsel(envdata,frames,0,1,envx(c));
            if max(curenv) > ymax
                ymax = max(curenv);
            end
        else           % minimum envelope
            curenv = matsel(envdata,frames,0,2,envx(c));
            if min(curenv) < ymin
                ymin = min(curenv);
            end
        end
        p=plot(times,curenv);

        % First, the outermost envelope
        if c == 1
            set(p, 'Color', colors{16,1},'LineWidth',2, 'LineStyle', '-');
        % After the second, each cluster
        elseif isempty(intersect(colorlog(:,1), envx(c)))
            if colorlog(1,1) == 0 && colorlog(1,2) == 0 % initial state
                colorlog(1,:) = [envx(c) 1];
            else
                colorlog(end+1,:) = [envx(c) size(colorlog,1)+1];
            end

            if     colorlog(end,2) < 13
                set(p, 'Color' ,linecolors{colorlog(end,2)},   'LineWidth',1, 'LineStyle', '-');
            elseif colorlog(end,2) < 25
                set(p, 'Color' ,linecolors{colorlog(end,2)-12},'LineWidth',3, 'LineStyle', ':');
            elseif colorlog(end,2) < 37
                set(p, 'Color' ,linecolors{colorlog(end,2)-24},'LineWidth',1, 'LineStyle', '--');
            elseif colorlog(end,2) < 49
                set(p, 'Color' ,linecolors{colorlog(end,2)-36},'LineWidth',2, 'LineStyle', '-.');
            elseif colorlog(end,2) < 61
                set(p, 'Color' ,linecolors{colorlog(end,2)-48},'LineWidth',2, 'LineStyle', '-');
            else   % this is for cluster numbers more than 60 
                set(p, 'Color', colors{16,1},                  'LineWidth',1, 'LineStyle', ':');
            end
        else
            [d, IA, IB] = intersect(colorlog(:,1), envx(c));
            if     IA < 13
                set(p, 'Color' ,linecolors{IA},   'LineWidth',1, 'LineStyle', '-');
            elseif IA < 25
                set(p, 'Color' ,linecolors{IA-12},'LineWidth',3, 'LineStyle', ':');
            elseif IA < 37
                set(p, 'Color' ,linecolors{IA-24},'LineWidth',1, 'LineStyle', '--');
            elseif IA < 49
                set(p, 'Color' ,linecolors{IA-36},'LineWidth',2, 'LineStyle', '-.');
            elseif IA < 61
                set(p, 'Color' ,linecolors{IA-48},'LineWidth',2, 'LineStyle', '-');
            else   % this is for cluster numbers more than 60 
                set(p, 'Color', colors{16,1}   ,'LineWidth',1, 'LineStyle', ':');
            end
        end
    end
end
set(gca,'FontSize',12,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fill specified component %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(g.fillclust)
    if ismember(g.fillclust, g.clust_grandERP(compx))
        [dummy c] = ismember(g.fillclust, g.clust_grandERP(compx));
        fprintf('filling the envelope of component %d\n',g.fillclust);
        mins = matsel(envdata,frames,0,2,compx(c)+1);
        p=fill([times times(frames:-1:1)],...
        [matsel(envdata,frames,0,1,compx(c)+1) mins(frames:-1:1)],[1 0 0]);
        % Overplot the data envlope again so it is not covered by the filled component
        p=plot(times,matsel(envdata,frames,0,1,1), 'k', 'LineWidth',2);% plot the max
        p=plot(times,matsel(envdata,frames,0,2,1), 'k', 'LineWidth',2);% plot the min
    else
        fprintf('cluster %d is not on the list for plotting\n',g.fillclust);
    end
end
clear clustlist clustn c dummy
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add vertical lines %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  

% plot vertical line at time zero
if g.timerange(1) <= 0 & g.timerange(2) >= 0
    vl=plot([0 0], [-1e10 1e10],'k');
    set(vl,'linewidth',2);
end

% if specified by option, plot specified vertical lines
if ~isempty(g.vert)
    for v=1:length(g.vert)
        vl=plot([g.vert(v) g.vert(v)], [-1e10 1e10],'k--');
        if any(g.limcontrib ~= 0) & v>= length(g.vert)-1;
            set(vl,'linewidth',LIMCONTRIBWEIGHT);
            set(vl,'linestyle',':');
        else
            set(vl,'linewidth',VERTWEIGHT);
            set(vl,'linestyle','--');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%%%%%%%%%%%%%%%%%%%% Extend y limits by 5% %%%%%%%%%%%%%%%%%%%%%%%%%%
%

datarange = ymax-ymin;
ymin = ymin-0.05*datarange;
ymax = ymax+0.05*datarange;
    
axis([pmin pmax ymin ymax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show pvaf/pv/rv/fillclust %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(g.sortvar,'pvaf')| strcmpi(g.sortvar,'pv')
    t = text(double(xmin+0.05*(xmax-xmin)), ...
        double(ymin+0.12*(ymax-ymin)), ...
        ['pvaf ' num2str(sumpvaf,'%4.2f') '%']);
    set(t,'fontsize',13,'fontweight','bold')
elseif strcmpi(g.sortvar,'rv') 
    t = text(double(xmin+0.05*(xmax-xmin)), ...
        double(ymin+0.12*(ymax-ymin)), ...
        ['rv ' num2str(sumpvaf,'%4.2f') '%']);
    set(t,'fontsize',13,'fontweight','bold')
end

if ~isempty(g.fillclust) & ismember(g.fillclust-1, envx)
    t = text(double(xmin+0.05*(xmax-xmin)), ...
        double(ymin+0.04*(ymax-ymin)), ...
        ['cluster ' num2str(g.fillclust) ' filled with red']);
    set(t,'fontsize',12,'fontweight','bold')
end

%
%%%%%%%%%%%%%%%%%%%%%% Label axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
set(axe,'Color',axcolor);
if strcmpi(g.xlabel, 'on')
    if ~xunitframes
        l= xlabel('Time (s)');
    else % xunitframes == 1
        l= xlabel('Data (time points)');
    end
    set(l,'FontSize',14,'FontWeight','Bold');
end
if strcmpi(g.ylabel, 'on')
    if strcmpi(g.envmode, 'avg')
        l=ylabel('Potential (uV)');
    else
        l=ylabel('RMS of uV');
    end;
    set(l,'FontSize',14,'FontWeight','Bold');
end
%
%%%%%%%%%%%%%% Draw maps and oblique/vertical lines %%%%%%%%%%%%%%%%%%%%%
%
% axall = axes('Units','Normalized','Position',pos,...
axall = axes('Position',pos,...
    'Visible','Off','Fontsize',16); % whole-figure invisible axes
axes(axall)
set(axall,'Color',axcolor);
axis([0 1 0 1])

width  = xmax-xmin;
pwidth  = pmax-pmin;
height = ymax-ymin;

if strcmpi(g.dispmaps, 'on')
    for t=1:ntopos % draw oblique lines from max env vals (or plot top) to map bases, in left to right order
        %
        %%%%%%%%%%%%%%%%%%% draw oblique lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if BOLD_COLORS==1
            linestyles = 1:ntopos;
        else
            linestyles = maporder;
        end
        axes(axall)
        axis([0 1 0 1]);
        set(axall,'Visible','off');
        maxenv = matsel(envdata,frames,plotframes(t),1,compx(t)+1);
        % max env val
        data_y = 0.6*(g.voffsets(t)+maxenv-ymin)/height;
        if (data_y > pos(2)+0.6*pos(4))
            data_y = pos(2)+0.6*pos(4);
        end
        
        % plot the oblique line
        l1 = plot([(plottimes(t)-pmin)/pwidth topoleft + 1/pos(3)*(t-1)*1.2*topowidth + (topowidth*0.6)], [data_y 0.68]); % 0.68 is bottom of topo maps
        
        % match cluster number and color using colorlog
        [d IA IB] = intersect(colorlog(:,1), compx(t)+1);
        if     IA < 13
            set(l1, 'Color' ,linecolors{IA},'LineWidth',1, 'LineStyle', '-');
        elseif IA < 25
            set(l1, 'Color' ,linecolors{IA-12},'LineWidth',3, 'LineStyle', ':');
        elseif IA < 37
            set(l1, 'Color' ,linecolors{IA-24},'LineWidth',1, 'LineStyle', '--');
        elseif IA < 49
            set(l1, 'Color' ,linecolors{IA-36},'LineWidth',2, 'LineStyle', '-.');
        elseif IA < 61
            set(l1, 'Color' ,linecolors{IA-48},'LineWidth',2, 'LineStyle', '-');
        else   % this is for cluster numbers more than 60 
            set(l1, 'Color', colors{16,1}   ,'LineWidth',1, 'LineStyle', ':');
        end
        
        hold on
        %
        %%%%%%%%%%%%%%%%%%%% add specified vertical lines %%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if g.voffsets(t) > 0
            l2 = plot([(plottimes(t)-xmin)/width (plottimes(t)-xmin)/width],[0.6*(maxenv-ymin)/height 0.6*(g.voffsets(t)+maxenv-ymin)/height]);
            if isempty(intersect(colorlog(:,1), compx(t)+1))
                if     t > 60
                    set(l2, 'Color', colors{16,1}         ,'LineWidth',1, 'LineStyle', ':');
                elseif t > 48
                    set(l2, 'Color' ,linecolors{mod(t,12)},'LineWidth',2, 'LineStyle', '-');
                elseif t > 36
                    set(l2, 'Color' ,linecolors{mod(t,12)},'LineWidth',2, 'LineStyle', '-.');
                elseif t > 24
                    set(l2, 'Color' ,linecolors{mod(t,12)},'LineWidth',1, 'LineStyle', '--');
                elseif t > 12
                    set(l2, 'Color' ,linecolors{mod(t,12)},'LineWidth',3, 'LineStyle', ':');
                elseif t > 0
                    set(l2, 'Color' ,linecolors{mod(t,12)} ,'LineWidth',1, 'LineStyle', '-');
                end
            else
                d = intersect(colorlog(:,1), compx(t)+1);
                if     d > 60
                    set(l2, 'Color' ,linecolors{colorlog(IA, 2)},'LineWidth',1, 'LineStyle', ':');
                elseif d > 48
                    set(l2, 'Color' ,linecolors{colorlog(IA, 2)},'LineWidth',2, 'LineStyle', '-');
                elseif d > 36
                    set(l2, 'Color' ,linecolors{colorlog(IA, 2)},'LineWidth',2, 'LineStyle', '-.');
                elseif d > 24
                    set(l2, 'Color' ,linecolors{colorlog(IA, 2)},'LineWidth',1, 'LineStyle', '--');
                elseif d > 12
                    set(l2, 'Color' ,linecolors{colorlog(IA, 2)},'LineWidth',3, 'LineStyle', ':');
                elseif d > 0
                    set(l2, 'Color' ,linecolors{colorlog(IA, 2)},'LineWidth',1, 'LineStyle', '-');
                end
            end
        end
        set(gca,'Visible','off');
        axis([0 1 0 1]);
    end % t
end; % if g.dispmaps == on

%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%
%

if strcmpi(g.dispmaps, 'on')

    % common scale for colors
    % -----------------------
    if strcmpi(g.actscale, 'on')
        maxvolt = 0;
        for n=1:ntopos
            maxvolt = max(max(abs(maxproj(:,n))), maxvolt);
        end;
    end;

    for t=1:ntopos % left to right order  (maporder)
        axt = axes('Units','Normalized','Position',...
            [pos(3)*topoleft+pos(1)+(t-1)*head_sep*topowidth pos(2)+0.66*pos(4) ...
            topowidth topowidth*head_sep]);
        axes(axt)                             % topoplot axes
        cla

        if g.gridind ~= 0
            tmp = zeros(67,67);
            tmp(:)=nan ;
            tmp(g.gridind) = maxproj(:,t);
            %tmp(g.gridind) = projERP{compx(t)}(:,plotframes(t));
        end
        
        figure(myfig);
        if strcmp(g.topotype, 'inst')
            toporeplot(tmp, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off');
        else % which is g.topytype == 'ave'
            toporeplot(STUDY.cluster(1, g.clust_grandERP(compx(t))).topo, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off');
        end

        axis square
        %
        %%%%%%%%%%%%% Scale colors %%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if strcmpi(g.actscale, 'on')
            caxis([-maxvolt maxvolt]);
        end;
        %
        %%%%%%%%%%%%%%%%%%%%%%%% label components %%%%%%%%%%%%%%%%%%%%%%%
        %
        if t==1
            chid = fopen('envtopo.labels','r');
            if chid <3,
                numlabels = 1;
            else
                fprintf('Will label scalp maps with labels from file %s\n','envtopo.labels');
                compnames = fscanf(chid,'%s',[4 MAXPLOTDATACHANS]);
                compnames = compnames';
                [r c] = size(compnames);
                for i=1:r
                    for j=1:c
                        if compnames(i,j)=='.',
                            compnames(i,j)=' ';
                        end;
                    end;
                end;
                numlabels=0;
            end
        end
        if numlabels == 1
            if ~isempty(complabels)
                complabel = complabels(t);
            else
                complabel = int2str(maporder(t));        % label comp. numbers
            end
        else
            complabel = compnames(t,:);              % use labels in file
        end
        text(0.00,0.80,complabel,'FontSize',14,...
            'FontWeight','Bold','HorizontalAlignment','Center');
        if strcmpi(g.sortvar, 'pvaf') | strcmpi(g.sortvar, 'pv')
            text(-0.6, -0.6, ['pvaf: ' sprintf('%6.2f', pvaf(t)) ] );
        elseif strcmpi(g.sortvar, 'rv')
            text(-0.4, -0.7, ['rv: ' sprintf('%6.2f', pvaf(t)) ] );
        end;
        % axt = axes('Units','Normalized','Position',[0 0 1 1],...
        axt = axes('Position',[0 0 1 1],...
            'Visible','Off','Fontsize',16);
        set(axt,'Color',axcolor);           % topoplot axes
        drawnow
    end

    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot a colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % axt = axes('Units','Normalized','Position',[.88 .58 .03 .10]);
    axt = axes('Position',[pos(1)+pos(3)*1.015 pos(2)+0.6055*pos(4) pos(3)*.02 pos(4)*0.09]);
    if strcmpi(g.actscale, 'on')
        h=cbar(axt, [1:64],[-maxvolt maxvolt],3);
    else
        h=cbar(axt);                        % colorbar axes
        set(h,'Ytick',[]);

        axes(axall)
        set(axall,'Color',axcolor);
        tmp = text(0.50,1.05,g.title,'FontSize',16,'HorizontalAlignment','Center','FontWeight','Bold');
        set(tmp, 'interpreter', 'none');
        text(1,0.68,'+','FontSize',16,'HorizontalAlignment','Center');
        % text(1,0.637,'0','FontSize',12,'HorizontalAlignment','Center','verticalalignment','middle');
        text(1,0.61,'-','FontSize',16,'HorizontalAlignment','Center');
    end;
    axes(axall)
    set(axall,'layer','top'); % bring component lines to top

end;
%
%%%%%%%%%%%%%%%%%%%%%%%%% turn on axcopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end;');

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function envdata = envelope(data, envmode)  % also in release as env()
if nargin < 2
    envmode = 'avg';
end;
if strcmpi(envmode, 'rms');
    warning off;
    negflag = (data < 0);
    dataneg = negflag.* data;
    dataneg = -sqrt(sum(dataneg.*dataneg,1) ./ sum(negflag,1));
    posflag = (data > 0);
    datapos = posflag.* data;
    datapos = sqrt(sum(datapos.*datapos,1) ./ sum(posflag,1));
    envdata = [datapos;dataneg];
    warning on;
else
    if size(data,1)>1
        maxdata = max(data); % max at each time point
        mindata = min(data); % min at each time point
        envdata = [maxdata;mindata];
    else
        maxdata = max([data;data]); % max at each time point
        mindata = min([data;data]); % min at each time point
        envdata = [maxdata;mindata];
    end
end;

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [STUDY, ALLEEG] = group_topo_now(STUDY, ALLEEG);
% this function is to separate topoall into number of groups
% written by Makoto Miyakoshi

% Is group the first or second variable?
if strcmp(STUDY.design(STUDY.currentdesign).variable(1,1).label, 'group')
    variable_group = 1;
else
    variable_group = 2;
end

% Obtain forward-reference; from STUDY.design(STUDY.currentdesign).cell to sets
% This is also an experiment to extract values from cell without using for_loop.
sets1 = cell2mat({STUDY.design(1,STUDY.currentdesign).cell.dataset}');
sets2 = strcat(STUDY.design(1,STUDY.currentdesign).cell.value);
sets2 = (sets2{1, variable_group})';
sets2 = str2num(sets2);
sets = cat(2, sets1, sets2); % column 1 for dataset numbers, column 2 for group numbers
clear sets1 sets2

% It seems that if values contained by cells embedded in cells, it is hard to avoid
% using for_loop (thus it is desireble to avoid such data structure).

% loop for all clusters
for cls = 2:length(STUDY.cluster) 
    if ~isempty(STUDY.cluster(1,cls).topoall) % this is for empty cluster 04/29/12 makoto
        % pick up set indices from the nth cluster
        ncomps = STUDY.cluster(1,cls).sets(1,:);

        for columnwise = 1:size(STUDY.cluster(1,cls).setinds, 1)
            for rowwise = 1:size(STUDY.cluster(1,cls).setinds, 2)
                for nsets = 1:length(STUDY.cluster(1,cls).setinds{columnwise,rowwise})
                    sets_grouped{columnwise, rowwise}(1, nsets) = STUDY.design(1,STUDY.currentdesign).cell(1, STUDY.cluster(1,cls).setinds{columnwise, rowwise}(1, nsets)).dataset;
                end
            end
        end

        % rotate if variable_group == 1;
        if variable_group == 1
            sets_grouped = sets_grouped';
        end

        sets_grouped_mat = cell2mat(sets_grouped);

        % reverse reference: from setinds to sets_divided_mat
        [a revref_group] = sort(sets_grouped_mat, 2);
        revref_group = revref_group(1,:);

        % create topoall_reordered
        for n = 1:size(STUDY.cluster(1,cls).sets, 2)
            topoall_reordered{1,n} = STUDY.cluster(1,cls).topoall{1,revref_group(n)};
        end

        % build topoall_group
        if variable_group == 1;
            ngroup = size(STUDY.cluster(1,cls).setinds, 1);
            for n = 1:ngroup
                ntopos(n) = length(STUDY.cluster(1,cls).setinds{n,1});
            end
        else
            ngroup = size(STUDY.cluster(1,cls).setinds, 2);
            for n = 1:ngroup
                ntopos(n) = length(STUDY.cluster(1,cls).setinds{1,n});
            end
        end

        n = 0;
        while n < ngroup
            n = n+1;
            if n == 1;
                topoall_group{1, n} = topoall_reordered(1:ntopos(n));
            else
                topoall_group{1, n} = topoall_reordered(ntopos(n-1)+1: ntopos(n-1)+ntopos(n));
            end
        end

        % locate the topoall_group to each of STUDY.cluster
        STUDY.cluster(1,cls).topoall_group = topoall_group;
        clear ncomps ntopos sets_* topoall_*
    end
end
