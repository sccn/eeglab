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
%
%   STUDY        = an EEGLAB STUDY structure containing EEG structures
%
%   ALLEEG       = the ALLEEG data structure; can also be an EEG dataset structure.
%
% Optional inputs:
%
%  'amplimits'   = [minuV maxuV]. {default: use data uV limits}
%
%  'baseline'    = [minms maxms] - a new baseline to remove from the grand
%                   and cluster ERPs.
%
%  'clustnums'   = [integer array] vector of cluster numbers to plot.  Else if
%                   int < 0, the number of largest contributing clusters to plot
%                   {default|[] -> 7}
%
%  'conditions'  = [integer array] vector of condition indices to plot
%
%  'diff'        = [condition1 group 1 condition 2 group 2] Perform subtractiono
%                   between conbination of condition/group. Ex. 'diff', [1 1 2 1]
%                   performs condition1_group1 minus condition2_group1.
%                   Must be provided with 4 numbers. If no condition/group, 1.
%
%  'fillclust'   = [integer] fill the numbered cluster envelope with red. {default|[]|0 -> no fill}
%
%  'fillcolor'   = [a b c] where a, b, c are =>0 and =<1. to create a color
%                   to fill the summed selected clusters. {dafault [0.875 0.875 0.875]} 
%
%  'limcontrib'  = [minms maxms]  time range (in ms) in which to rank cluster contributions
%                   (boundaries = thin dotted lines) {default|[]|[0 0] -> plotting limits}
%
%  'onlyclus'    = [ 'on' | 'off'] dataset components to include in the grand ERP.
%                  'on' will include only the components that were part of the
%                   clustering. For example, if components were rejected from
%                   clustering because of high dipole model residual variance,
%                   don't include their data in the grand ERP.
%                  'off' will include all components in the datasets except
%                   those in the subtructed ('subclus') clusters {default 'on'}.
%
%  'sortvar'     = ['maxp'|'pvaf'|'ppaf'|'rltp'|'area']
%                  'maxp', maximum power
%                    maxp(comp) = max(sum(cluster(:,selectedtimewindow)^2));
%                  'pvaf', sort components by percent variance accounted for (eeg_pvaf())
%                    pvaf(comp) = 100-100*mean(var(data - back_proj))/mean(var(data));
%                  'ppaf', sort components by percent power accounted for (ppaf) 
%                    ppaf(comp) = 100-100*Mean((data - back_proj).^2)/Mean(data.^2);
%                  'rltp', sort components by relative power 
%                    rltp(comp) = 100*Mean(back_proj.^2)/Mean(data.^2);
%                  'area', sort components by enveloped area ratio
%                    area(comp) = 100*(area of envelope by a selected cluster)/(area of outermost envelope)
%
%  'sortvarnorm' = ['on'|'off'] is normalizes sortvar so that the sum of sortvar is 100%  
%
%  'subclus'     = [integer array] vector of cluster numbers to omit when  computing
%                    the ERP envelope of the data (e.g., artifact
%                    clusters). By default, these clusters are excluded from
%                   the grandERP envelope. See the next option 'onlyclus'.
%
%  'sumenvfill'  = ['selective'|'all'|'off'] fill or not the envelope of the
%                    summed selected cluster. 'select' fills only limcontrib
%                    time range. 'all' fills all timerange. 'off' does not fill.
%                    See also 'fillcolor' so choose a filling color. {default: 'select'}
%
%  'timerange'   = data epoch start and end input latencies (in ms)
%                   {default: from 'limits' if any}
%
%  'topotype'    = ['inst'|'ave'] If 'inst', show the instantaneous map at the specific timepoint specified by the line.
%                   If 'ave', show that of the mean which is the same as the
%                   stored topomap. {default: 'inst'}
%
%  'vert'        = vector of times (in ms) at which to plot vertical dashed lines
%                   {default|[] -> none}
%
% See also: eegplugin_std_envtopo, pop_std_envtopo, std_envtopo, envtopo

% Author: Makoto Miyakoshi, Hilit Serby, Arnold Delorme, Scott Makeig
% History:
% 01/30/2013 ver 6.1 by Makoto. group_topo_now related bug fixed.
% 01/28/2013 ver 6.0 by Makoto. Combined groups supported. group_topo_now simplified.
% 10/10/2012 ver 5.3 by Makoto. 'sortvar' 'sortvarnorm' 'fillcolor' 'sumenvfill' added.
% 04/29/2012 ver 5.2 by Makoto. Bug fixed. STUDY.design(STUDY.currentdesign)
% 04/24/2012 ver 5.1 by Makoto. Revived 'amplimit' 'fillclust' options, and vertical doted lines for limcontrib range. Default changed into 'pvaf' from 'rv'. 
% 02/28/2012 ver 5.0 by Makoto. Values are now normalized by the number of subjects participating to the cluster. Group difference is better addressed when performing subtraction across groups. A group ratio is output too.
% 02/24/2012 ver 4.5 by Makoto. New color schema to improve color readability & shared line colors across multiple plots
% 02/22/2012 ver 4.4 by Makoto. Now plot title can accept numbers
% 01/24/2012 ver 4.3 by Makoto. Output added for the function.
% 01/17/2012 ver 4.2 by Makoto. Improved diff option; now mutually exclusive with normal plot. Diff title improved too. 
% 01/06/2012 ver 4.1 by Makoto. Time range options fixed. STUDY.design selection fixed. Wrong labels fixed. (Thanks to Agatha Lenartowicz!) 'sortvar' and 'topotype' implemented in the main function. Help edited.
% 12/30/2011 ver 4.0 by Makoto. Main function 100% redesigned. Full support for Study design.
% 12/26/2011 ver 3.2 by Makoto. Formats corrected.
% 12/25/2011 ver 3.1 by Makoto. Bug fixed about choosing wrong clusters under certain conditions. 
% 10/12/2011 ver 3.0 by Makoto. Completed the alpha version. 
% 10/11/2011 ver 2.2 by Makoto. Time range, group with no entry fixed
% 10/10/2011 ver 2.1 by Makoto. Two types of scalp topos can be presented.
% 10/08/2011 ver 2.0 by Makoto. Multiple groups supported (but one group one time using option). Result topos are now retrieved from the stored ones in STUDY 
% 08/01/2011 ver 1.0 by Makoto. Updated the original script to read data from memory.

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
              'sortvar'       'string'   {'maxp', 'pvaf', 'ppaf', 'rltp', 'area'}         'pvaf';...
              'sortvarnorm'   'string'   {'on', 'off'}                                     'off';...
              'sumenvfill'    'string'   {'selective', 'all', 'off'}                 'selective';...
              'topotype'      'string'   {'inst', 'ave'}                                  'inst';...
              'fillcolor'     'real'     []                                [0.875  0.875  0.875];...
              'fillclust'     'integer'  []                                                   []});
          
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
    if isfield(STUDY.cluster, 'topoall_group')
        for n = 2:length(STUDY.cluster)
            if isempty(STUDY.cluster(n).topoall_group)
                [STUDY, ALLEEG] = group_topo_now(STUDY, ALLEEG);
                break
            end
        end
    else
        [STUDY, ALLEEG] = group_topo_now(STUDY, ALLEEG);
    end
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
            % Check if conditions/group combined
            tmpCond  = STUDY.design(STUDY.currentdesign).variable(1,1).value{1, columnwise};
            tmpGroup = STUDY.design(STUDY.currentdesign).variable(1,2).value{1, rowwise};
            if iscell(tmpCond);  tmpCond  = cell2mat(tmpCond);  end
            if iscell(tmpGroup); tmpGroup = cell2mat(tmpGroup); end
            arglist.title = [num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(tmpCond) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(tmpGroup)];
            colorlog = envtopo_plot(STUDY, ALLEEG, squeeze(outermostenv(:,arglist.timepoints,1,columnwise,rowwise)), squeeze(topoerpconv_pile(:,arglist.timepoints,:,columnwise,rowwise)), colorlog, arglist);
        end
    end
end
clear columnwise rowwise tmpCond tmpGroup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Else, plot difference envtopo    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if ~isempty(arglist.diff)
    
    % First, plot title is determined
    tmpCond1  = STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.diff(1,1)};
    tmpCond2  = STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.diff(1,3)};
    tmpGroup1 = STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.diff(1,2)};
    tmpGroup2 = STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.diff(1,4)};
    if iscell(tmpCond1);  tmpCond1  = cell2mat(tmpCond1);  end
    if iscell(tmpCond2);  tmpCond2  = cell2mat(tmpCond2);  end        
    if iscell(tmpGroup1); tmpGroup1 = cell2mat(tmpGroup1); end
    if iscell(tmpGroup2); tmpGroup2 = cell2mat(tmpGroup2); end
    
    % Group only
    if     length(STUDY.design(STUDY.currentdesign).variable(1,1).value) == 1 && arglist.diff(1) == 1 && arglist.diff(3) == 1
        arglist.title = [num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(tmpGroup1) ' minus ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(tmpGroup2)];
    % Condition only
    elseif length(STUDY.group) == 1 && arglist.diff(2) == 1 && arglist.diff(4) == 1 % Condition only if no group
        arglist.title = [num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(tmpCond1)  ' minus ' num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(tmpCond2)];
    % if both condition and group
    else   
        arglist.title = [num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(tmpCond1) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(tmpGroup1) ' minus ' num2str(STUDY.design(STUDY.currentdesign).variable(1,1).label) ' ' num2str(tmpCond2) ' ' num2str(STUDY.design(STUDY.currentdesign).variable(1,2).label) ' ' num2str(tmpGroup2)];
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





%
% envtopo_plot() - Plot the envelope of a data epoch, plus envelopes and scalp maps of specified
%             or largest-contributing components. If a 3-D input matrix, operates on the
%             mean of the data epochs. Click on individual axes to examine them in detail.
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

if ~xunitframes
    fprintf('  in the interval %3.0f ms to %3.0f ms.\n',1000*times(limframe1),1000*times(limframe2));
end

% sort components by maximum mean back-projected power 
% in the 'limcontrib' time range: mp(comp) = max(Mean(back_proj.^2));
% where back_proj = comp_map * comp_activation(t) for t in 'limcontrib'

%
%%%%%% Calculate sortvar, used to sort the components %%%%%%%%%%%
%
for c = 1:length(projERP)
    if     strcmpi(g.sortvar,'maxp')  % Maximum Power of backproj
        sortvar(c) = max(mean(projERP{c}(:,limframe1:limframe2).*projERP{c}(:,limframe1:limframe2)));
        fprintf('%s maximum mean power of back-projection: %g\n',g.clustlabels{c},sortvar(c));

    elseif strcmpi(g.sortvar,'pvaf')   % Percent Variance
        vardat = var(reshape(grandERP(:,limframe1:limframe2),1,nvals));
        difdat = grandERP(:,limframe1:limframe2)-projERP{c}(:,limframe1:limframe2);
        difdat = reshape(difdat,1,nvals);
        sortvar(c) = 100-100*(var(difdat)/vardat); %var of diff div by var of full data
        fprintf('%s percent variance accounted for(pvaf): %2.2f%%\n',g.clustlabels{c},sortvar(c));
        
    elseif strcmpi(g.sortvar,'ppaf')    % Percent Power
        powdat = mean(mean(grandERP(:,limframe1:limframe2).^2));
        sortvar(c) = 100-100*mean(mean((grandERP(:,limframe1:limframe2)...
                    -projERP{c}(:,limframe1:limframe2)).^2))/powdat;
        fprintf('%s percent power accounted for(ppaf): %2.2f%%\n',g.clustlabels{c},sortvar(c));

    elseif strcmpi(g.sortvar,'rltp')    % Relative Power
        powdat = mean(mean(grandERP(:,limframe1:limframe2).^2));
        sortvar(c) = 100*mean(mean((projERP{c}(:,limframe1:limframe2)).^2))/powdat;
        fprintf('%s relative power of back-projection: %2.2f%%\n',g.clustlabels{c},sortvar(c));
        
    elseif strcmpi(g.sortvar,'area')    % Enveloped Area ratio
        outermostenv      = envelope(grandERP(:,limframe1:limframe2), g.envmode);
        outermostenv_area = sum(abs(outermostenv(1,:)))+sum(abs(outermostenv(2,:)));
        tmp_projERP       = envelope(projERP{c}(:,limframe1:limframe2), g.envmode);
        tmp_projERP_area  = sum(abs(tmp_projERP(1,:)))+sum(abs(tmp_projERP(2,:)));
        sortvar(c)        = 100*tmp_projERP_area/outermostenv_area;
        fprintf('%s enveloped area ratio : %2.2f%%\n',g.clustlabels{c},sortvar(c));
    end
end

%%%%%% Normalize the sortvar to 100 if requested %%%%%%
if strmatch(g.sortvarnorm, 'on')
    sortvar = sortvar*1/sum(sortvar)*100;
    fprintf('\n')
    fprintf('sortvarnorm is on: sortvar is normalized to 100%% in the plot \n');
end

%
%%%%%%%%%%%%%%%%%%%%%%%%% Sort by max variance in data %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[sortvar,compx] = sort(sortvar);  % sort clustnums on max sortvar
sortvar = sortvar(ncomps:-1:1);  % reverse order of sort (max:min)       
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
sortvar = sortvar(1:ntopos);
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
sortvar = sortvar(ifx);
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

fprintf('    Component sortvar in interval:  ');
for t=1:ntopos
    fprintf('%s ',num2str(sortvar(t)));
end
fprintf('\n');

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

sumproj = zeros(size(projERP{1}));
for n = 1:ntopos
    sumproj = sumproj + projERP{compx(n)}; % add up all cluster projections
end

% calculate summed sortvar of selected clusters
if     strcmpi(g.sortvar,'maxp')  % Maximum Power of backproj
    selectclusvar = max(mean(sumproj(:,limframe1:limframe2).*sumproj(:,limframe1:limframe2)));

elseif strcmpi(g.sortvar,'pvaf')   % Percent Variance
    vardat = var(reshape(grandERP(:,limframe1:limframe2),1,nvals));
    difdat = grandERP(:,limframe1:limframe2)-sumproj(:,limframe1:limframe2);
    difdat = reshape(difdat,1,nvals);
    selectclusvar = 100-100*(var(difdat)/vardat); %var of diff div by var of full data

elseif strcmpi(g.sortvar,'ppaf')    % Percent Power
    powdat = mean(mean(grandERP(:,limframe1:limframe2).^2));
    selectclusvar = 100-100*mean(mean((grandERP(:,limframe1:limframe2)...
                -sumproj(:,limframe1:limframe2)).^2))/powdat;

elseif strcmpi(g.sortvar,'rltp')    % Relative Power
    powdat = mean(mean(grandERP(:,limframe1:limframe2).^2));
    selectclusvar = 100*mean(mean((sumproj(:,limframe1:limframe2)).^2))/powdat;

elseif strcmpi(g.sortvar,'area')    % Enveloped Area ratio
    outermostenv      = envelope(grandERP(:,limframe1:limframe2), g.envmode);
    outermostenv_area = sum(abs(outermostenv(1,:)))+sum(abs(outermostenv(2,:)));
    tmp_projERP       = envelope(sumproj(:,limframe1:limframe2), g.envmode);
    tmp_projERP_area  = sum(abs(tmp_projERP(1,:)))+sum(abs(tmp_projERP(2,:)));
    selectclusvar        = 100*tmp_projERP_area/outermostenv_area;
end

fprintf('Summed cluster sortvar in interval [%s %s] ms: %4.2f%',...
    int2str(1000*times(limframe1)),int2str(1000*times(limframe2)), selectclusvar);
fprintf('\n')



if isempty(g.amplimits)
    envx = [1,compx+1];
    ymin = min(matsel(envdata,frames,0,1,envx(1)));
    ymax = max(matsel(envdata,frames,0,2,envx(1)));
else
    ymin = g.amplimits(1);
    ymax = g.amplimits(2);
end

%
% Plot the summed projection filled
%
sumenv = envelope(sumproj, g.envmode);
if     strcmpi(g.sumenvfill,'selective')
    filltimes = times(limframe1:limframe2);
    mins = matsel(sumenv,frames,0,2,0);
    p=fill([filltimes filltimes(end:-1:1)],...
         [matsel(sumenv,frames,limframe1:limframe2,1,0) mins(limframe2:-1:limframe1)], g.fillcolor);
    set(p,'EdgeColor', g.fillcolor);
    hold on
    % redraw outlines
    p=plot(times,matsel(envdata,frames,0,1,1),'Color', colors{16,1});% plot the max
    set(p,'LineWidth',2);
    p=plot(times,matsel(envdata,frames,0,2,1),'Color', colors{16,1});% plot the min
    set(p,'LineWidth',2);
    
elseif strcmpi(g.sumenvfill,'all')
    mins = matsel(sumenv,frames,0,2,0);
    p=fill([times times(frames:-1:1)],...
        [matsel(sumenv,frames,0,1,0) mins(frames:-1:1)], g.fillcolor);
    set(p,'EdgeColor', g.fillcolor);
    hold on
    % redraw outlines
    p=plot(times,matsel(envdata,frames,0,1,1),'Color', colors{16,1});% plot the max
    set(p,'LineWidth',2);
    p=plot(times,matsel(envdata,frames,0,2,1),'Color', colors{16,1});% plot the min
    set(p,'LineWidth',2);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show sortvar and fillclust %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = text(double(xmin+0.05*(xmax-xmin)), double(ymin+0.12*(ymax-ymin)), ...
    [g.sortvar ' ' num2str(selectclusvar,'%4.2f') '%']);
set(t,'fontsize',13,'fontweight','bold')

if ~isempty(g.fillclust) & ismember(g.fillclust-1, envx)
    t = text(double(xmin+0.05*(xmax-xmin)), ...
        double(ymin+0.04*(ymax-ymin)), ...
        ['cluster ' num2str(g.fillclust) ' filled with a color']);
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
        
        text(-0.6, -0.6, [g.sortvar ': ' sprintf('%6.2f', sortvar(t)) ] );
        
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



function [STUDY, ALLEEG] = group_topo_now(STUDY, ALLEEG)
disp('Loading and separating topographs (only once)...')
var2Len = size(STUDY.design(STUDY.currentdesign).variable(2).value, 2);
for cls = 2:length(STUDY.cluster)
    for var2 = 1:var2Len
        for icLen = 1:length(STUDY.cluster(1,cls).setinds{1,var2})
            tmpAllinds = STUDY.cluster(1,cls).allinds{1,var2}(icLen);
            tmpSetinds = STUDY.cluster(1,cls).setinds{1,var2}(icLen);
            tmpDataset = STUDY.design(STUDY.currentdesign).cell(tmpSetinds).dataset;
            STUDY.cluster(1,cls).topoall_group{1,var2}{1,icLen} = std_readtopo(ALLEEG, tmpDataset, tmpAllinds);
        end
    end
end