% std_envtopo() - Creates an envtopo() image for a STUDY set using component cluster
%                 contributions instead of individual components.  Plots the envelope
%                 of the data epoch grand mean ERP, plus envelopes and average scalp maps
%                 for specified or largest-contributing clusters for each condition.
%                 Click on individual axes to examine them in detail (using axcopy()).
%                 See envtopo() for further details.
% Usage:
%                >> std_envtopo(STUDY, ALLEEG, 'key1', 'val1', ...);
% Inputs:
%   STUDY        = an EEGLAB STUDY structure containing EEG structures
%   ALLEEG       = the ALLEEG data structure; can also be an EEG dataset structure.
%
% Optional inputs:
%  'clusters'    = [integer array] vector of cluster numbers. If one cluster, plots the
%                    cluster contribution to the data envelope. If multiple clusters,
%                    selects the largest contributing clusters from within the
%                  'limcontrib' region (see below) and plots the envelopes of
%                    their contributions {default| [] -> 'all'}
%  'conditions'  = [integer array] vector of condition indices to plot
%  'subclus'     = [integer array] vector of cluster numbers to omit when  computing
%                    the ERP envelope of the data (e.g., artifact clusters)
%                    {default|[] -> none}
%  'env_erp'     = ['contrib' | 'all']
%                  'contrib' - > If one cluster, the grand ERP envelope includes
%                   only the datasets that are part of that cluster.
%                  'all' -> Grand ERP envelope includes all datasets in STUDY.
%                   If multiple clusters, this is the only option possible.
%  'only_clust'  = [ 'on' | 'off'] dataset components to include in the grand ERP.
%                  'on' will include only the components that were part of the
%                   clustering. For example, if components were rejected from
%                   clustering because of high dipole model residual variance,
%                   don't include their data in the grand ERP.
%                  'off' will include all components in the datasets except
%                   those in the subtructed ('subclus') clusters {default 'off'}.
%  'baseline'    = [minms maxms] - a new baseline to remove from the grand
%                   ERP and cluster ERP contributions.
%  'diff'        = [condition1 condition2] the numbers of two conditions.
%                   Plots an additional figure with the difference of the two conditions.
%  'clustnums'   = [integer array] vector of cluster numbers to plot.  Else if
%                   int < 0, the number of largest contributing clusters to plot
%                   {default|[] -> 7}
%  'timerange'   = data epoch start and end input latencies (in ms)
%                   {default: from 'limits' if any}
%  'limits'      = 0 or [minms maxms] or [minms maxms minuV maxuV]. Specify start/end plot
%                   x-axis limits (in ms) and min/max y-axis limits (in uV). If 0, or if both
%                   minmx & maxms == 0 -> use latencies from 'timerange' (else, 0:frames-1).
%                   If both minuV and maxuV == 0 -> use data uV limits {default: 0}
%  'limcontrib'  = [minms maxms]  time range (in ms) in which to rank cluster contributions
%                   (boundaries = thin dotted lines) {default|[]|[0 0] -> plotting limits}
%  'vert'        = vector of times (in ms) at which to plot vertical dashed lines
%                   {default|[] -> none}
%
% See also: envtopo()

% $Log: not supported by cvs2svn $
% Revision 1.24  2008/01/17 21:54:36  elizabeth
% clusterp.erp -> clusterp.erpdata
%
% Revision 1.23  2007/09/11 10:59:41  arno
% same as 1.20
%
% Revision 1.22  2007/09/11 10:44:52  arno
% nothing
%
% Revision 1.20  2007/08/06 22:18:23  arno
% do not use std_clustread any more
%
% Revision 1.19  2007/08/05 20:28:15  arno
% debug conditions
%
% Revision 1.18  2007/08/05 20:27:09  arno
% debug conditions etc...
%
% Revision 1.17  2007/08/05 20:26:06  arno
% version with conditions
%
% Revision 1.15  2007/05/26 03:40:34  toby
% Corrected title in case of a single cluster, removed obsolete commented code
%
% Revision 1.14  2007/05/26 03:37:39  toby
% diff option now displays all cluster topoplots and title is meaningful
%
% Revision 1.13  2007/05/25 03:44:51  toby
% RCS commenting added
%

function std_envtopo(STUDY, ALLEEG, varargin)
icadefs;
if nargin < 2
    help std_envtopo;
    return
end

if mod(nargin,2) % if not an even number of arguments
    error('std_envtopo: Input argument list must be pairs of: ''keyx'', ''valx'' ');
end

conditions   = [];       % default: all conditions
clusters     = [];       % default: use all clusters
subclus      = [];       % default: omit no clusters from the ERP envelope
e_options{1} = 'env_erp';
e_options{2} = 'contrib';% default: grand ERP includes only the components that
% contribute to the dataset; otherwise 'all' datasets in STUDY.
e_options{3} = 'only_clust'; % Include only those components that were part of...
% the pre-clustering data when computing the grand ERP
e_options{4} = 'off'; % Default is off: include all components in the datasets
% except those in the subtracted ('subclus') clusters
p_options = {};       % Default: other options off

for k = 1:2:(nargin-2)
    switch varargin{k}
        case 'clusters'
            clusters = varargin{k+1};
        case 'conditions'
            conditions = varargin{k+1};
        case 'env_erp'
            e_options{2} = varargin{k+1};
        case 'only_clust'
            e_options{4} = varargin{k+1};
        case 'subclus'
            subclus = varargin{k+1};
            e_options{end+1} = 'subclus';
            e_options{end+1} = subclus;
        case  'limcontrib'
            p_options{end+1} = 'limcontrib';
            p_options{end+1} = varargin{k+1};
        case  'limits'
            limits = varargin{k+1};
            if length(limits) == 4
                p_options{end+1} = 'limits';
                p_options{end+1} = limits;
            end
        case 'baseline'
            t = ALLEEG(STUDY.datasetinfo(STUDY.setind(1)).index).times;
            % convert baseline from ms to indexes
            maxind = max(find(t <= varargin{k+1}(end)));
            minind = min(find(t >= varargin{k+1}(1)));
            baseline = (minind:maxind);
            e_options{end+1} = 'baseline';
            e_options{end+1} = baseline;
        case 'vert'
            vert = varargin{k+1};
            p_options{end+1} = 'vert';
            p_options{end+1} = vert;
        case 'clustnums'
            p_options{end+1} = 'clustnums';
            p_options{end+1} = varargin{k+1};
        case 'diff'
            diffc = varargin{k+1};
        case 'timerange'
            timerange = varargin{k+1};
            p_options{end+1} = 'timerange';
            p_options{end+1} = timerange;
    end
end

if ~exist('timerange')
    timerange = [ALLEEG(STUDY.datasetinfo(STUDY.setind(1)).index).times([1 end])];
    p_options{end+1} = 'timerange';
    p_options{end+1} = timerange;
end

if isempty(clusters) | strcmpi(clusters,'all') % Default to all clusters in STUDY
    clusters = [];
    for k = 2:length(STUDY.cluster)
        if ~strncmpi('Notclust',STUDY.cluster(k).name,8)
            clusters = [clusters k];
        end
    end
end

% If some of the requested clusters-to-subtruct are in clusters, remove them.

clusters = setdiff(clusters,subclus);

%if length(clusters) > 1
%    e_options{2} = 'all';
%end

Ncond = length(STUDY.condition); % number of conditions
if Ncond == 0
    Ncond = 1;
end
if isempty(conditions)
    conditions = 1:Ncond;
end

for n = conditions
    %
    % Compute the grand mean ERP envelope of the cluster (for a specific condition).
    %
    fprintf('\n Computing grand ERP for condition %d.', n);
    [grandERP, set_len, STUDY, ALLEEG] = std_granderp(STUDY, ALLEEG, 'clusters', clusters, 'condition', n, e_options{:});
    grandERPtot{n} = grandERP;
    if ~exist('limits') | (length(limits)  == 2) % make limits same all conditions
        tmpmin = min(min(grandERP));
        tmpmax = max(max(grandERP));
        datarange = tmpmax-tmpmin;
        tmpmin = tmpmin-0.05*datarange;
        tmpmax = tmpmax+0.05*datarange;
        if n == conditions(1)
            ymin = tmpmin;
            ymax = tmpmax;
        else
            ymin = min(tmpmin,ymin);
            ymax = max(tmpmax,ymax);
        end
        if n == conditions(end)
            p_options{end+1} = 'limits';
            if ~exist('limits')
                p_options{end+1} = [timerange ymin ymax];
            else
                p_options{end+1} = [limits ymin ymax];
            end
        end
    end
    clear grandERP;
end

for n = conditions
    %
    % Compute the grand mean ERP envelope of the cluster
    % (for a specific condition).
    %
    for cls = 1:length(clusters)
        len = length(STUDY.cluster(clusters(cls)).comps);
        try
            % clustscalp = std_clustread(STUDY, ALLEEG, clusters(cls),'scalp'); % read scalp maps from cluster
            [tmp clustscalp] = std_readtopoclust(STUDY, ALLEEG, clusters(cls)); % read scalp maps from cluster
            clustscalp = clustscalp{1};
        catch,
            warndlg2([ 'Some topoplot information is missing, aborting'] , 'Abort - std_envtopo()' );
            return;
        end
        %clusterp = std_clustread(STUDY, ALLEEG, clusters(cls),'erp', n);

        %[tmp clusterp] = std_readdata(STUDY, ALLEEG, 'clusters', clusters(cls), 'infotype', 'erp');
        for k = 1:len
            dat  = STUDY.cluster(clusters(cls)).sets(n,k);
            comp = STUDY.cluster(clusters(cls)).comps(k);
            clusterp.erp{k} = std_readerp(ALLEEG, dat, comp);
            %clusterp.erp{k} = clusterp.erpdata{n}(:,k); % only works if 1 group
        end;
        if exist('baseline')
            for k = 1:length(clusterp.erp)
                clusterp.erp{k} = rmbase(clusterp.erp{k},...
                    ALLEEG(STUDY.datasetinfo(STUDY.setind(1)).index).pnts,baseline);
            end;
        end
        
        %
        % Compute grand mean back projection ERP for the cluster
        %
        projERP = 0;
        fprintf('\n Computing grand ERP projection of cluster %d: ', (clusters(cls)) );
        val_ind = find(~isnan(clustscalp.topo{1}(:))); % find non-NAN values
        for k = 1:len
            tmp = clustscalp.topo{k}(val_ind);
            projERP = projERP + tmp*clusterp.erp{k};
            fprintf('.');
        end
        tot_projERP{cls} = projERP/set_len;
        clust_names{cls} = [STUDY.cluster(clusters(cls)).name];
        numind = strfind(STUDY.cluster(clusters(cls)).name,' ')+1; % find the cluster number
        clus_ind{cls} = [STUDY.cluster(clusters(cls)).name(numind(end):end) ', '];
    end
    p_options{end+1} = 'clustlabels';
    p_options{end+1} = clust_names;
    p_options{end+1} = 'gridind';
    p_options{end+1} = val_ind;
    p_options{end+1} = 'sortvar';
    p_options{end+1} = 'pvaf';

    %
    % Compute the grand mean ERP envelope of the cluster (for a specific condition).
    %
    if cls == 1
        figure; set(gcf,'Color', BACKCOLOR);
        orient landscape;
        %
        % To save memory, once computation of a condition is complete
        % its grand ERP is cleared grandERPtot(1), which makes the current
        % condition grand ERP always in index 1
        %
        envtopo_plot(grandERPtot{1},tot_projERP,'envmode' ,'avg', 'fillcomp', 1, ...
            'dispmaps', 'on', 'title', ...
            ['Projected cluster: ' clus_ind{:} ' ' STUDY.condition{n}], p_options{:} );
    else
        figure; set(gcf,'Color', BACKCOLOR);
        orient landscape;
        envtopo_plot(grandERPtot{1},tot_projERP,'envmode' ,'avg', ...
            'dispmaps', 'on', 'title', [ STUDY.condition{n}], p_options{:} );
    end

    if exist('diffc')
        if  ~exist('diff_h')
            diff_h = figure;
            orient landscape;
            set(gcf,'Color', BACKCOLOR);
        end
        if n == diffc(1)
            figure(diff_h); subplot(3,1,1)
            envtopo_plot(grandERPtot{1},tot_projERP,'envmode' ,'avg', 'dispmaps', 'off', p_options{:},...
                'title', STUDY.condition{n} ,'xlabel', 'off' );
            erp1 = grandERPtot{1}; % The grand ERP
            proj1 = tot_projERP;
            clear  tot_projERP
        end
        if n == diffc(2)
            figure(diff_h); subplot(3,1,2)
            envtopo_plot(grandERPtot{1},tot_projERP,'envmode' ,'avg', 'dispmaps', 'off', p_options{:}, ...
                'title', STUDY.condition{n} ,'xlabel', 'off' );
            erp2 = grandERPtot{1}; % The grand ERP
            proj2 = tot_projERP;
            clear  tot_projERP
        end
        if exist('erp1') & exist('erp2')
            for cls = 1:length(clusters) 
                diff_proj{cls}  = proj1{cls} - proj2{cls};
            end
            clear proj1 proj2
            %plot difference
            figure(diff_h);
            subplot(3,1,3);
            envtopo_plot(erp1-erp2,diff_proj,'envmode' ,'avg', 'dispmaps','off',...
                'title','Difference',p_options{:} );
            
            figure
            if cls == 1
                envtopo_plot(erp1-erp2,diff_proj,'envmode' ,'avg', 'fillcomp', 1, ...
                    'title',['Difference Between Conditions ',int2str(diffc(1)),' and ',int2str(diffc(2)),'.'],...
                    'dispmaps', 'on', p_options{:} );
            else
                envtopo_plot(erp1-erp2,diff_proj,'envmode' ,'avg', 'dispmaps','on',...
                    'title',['Difference Between Conditions ',int2str(diffc(1)),' and ',int2str(diffc(2)),'.'],...
                    p_options{:} );
            end
            orient landscape;
            clear erp1 erp2  % remove variables so this loop runs only once
        end
    end
    grandERPtot(1) = [];
end

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
%  projERP          = A cell array of the projected ERPs of the desired components, each cell size is (chans,frames),
%                              (see std_granderp() and std_envtopo() for details).
%
% Optional inputs:
%  'clustnums'  = [integer array] vector of clusters numbers to plot {default|0 -> all}
%                  Else if int < 0, the number of largest contributing clusters to plot
%                  {default|[] -> 7}
%  'timerange' = start and end input data latencies (in ms) {default: from 'limits' if any}
%  'limits'    = 0 or [minms maxms] or [minms maxms minuV maxuV]. Specify start/end plot
%                  (x) limits (in ms) and min/max y-axis limits (in uV). If 0, or if both
%                  minmx & maxms == 0 -> use latencies from 'timerange' (else 0:frames-1).
%                  If both minuV and maxuV == 0 -> use data uV limits {default: 0}
%  'limcontrib' = [minms maxms]  time range (in ms) in which to rank component contribution
%                  (boundaries shown with thin dotted lines) {default|[]|[0 0] -> plotting limits}
%  'sortvar'    = ['mv'|'pv'|'rv'] if 'mv', sort components by back-projected variance; if 'pv',
%                  sort by percent variance accounted for (pvaf). If 'rv', sort by relative
%                  variance. Here:
%                   pvaf(component) = 100-100*variance(data-component))/variance(data)
%                   rv(component)   = 100*variance(component)/variance(data) {default: 'mv'}
%  'title'      = [string] plot title {default|[] -> none}
%  'plotchans'  = [integer array] data channels to use in computing contributions and envelopes,
%                  and also for making scalp topo plots {default|[] -> all}
%  'voffsets'   = [float array] vertical line extentions above the data max to disentangle
%                  plot lines (left->right heads, values in y-axis units) {def|[] -> none}
%  'colors'     = [string] filename of file containing colors for envelopes, 3 chars
%                  per line, (. = blank). First color should be "w.." (white)
%                  Else, 'bold' -> plot default colors in thick lines.
%                  {default|[] -> standard Matlab color order}
%  'fillcomp'   = int_vector>0 -> fill the numbered component envelope(s) with
%                  solid color. Ex: [1] or [1 5] {default|[]|0 -> no fill}
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

function [compvarorder,compvars,compframes,comptimes,compsplotted,pvaf] = envtopo_plot(grandERP,projERP,varargin)

pvaf = []; % note previous arguments 'on' -> 'pv', 'off' -> 'rv'
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

if nargin < 2
    help envtopo_plot
    return
end

if nargin <= 2 | isstr(varargin{1})
    % 'key' 'val' sequences
    fieldlist = { 'title'         'string'   []                       '';
        'clustlabels'    'cell'    {}                      {};
        'limits'        'real'     []                       0;
        'timerange'     'real'     []                       [];
        'voffsets'      'real'     []                       [] ;
        'vert'          'real'     []                       [] ;
        'fillcomp'      'integer'  []                       0 ;
        'colorfile'     'string'   []                       '' ;
        'colors'        'string'   []                       '' ;
        'clustnums'      'integer'  []                       -7;
        'dispmaps'      'string'   {'on' 'off'}             'on';
        'xlabel'      'string'   {'on' 'off'}             'on';
        'ylabel'      'string'   {'on' 'off'}             'on';
        'pvaf'          'string'   {'mv' 'rv' 'pv' 'pvaf' 'on' 'off' ''}  '';
        'envmode'       'string'   {'avg' 'rms'}            'avg';
        'sortvar'       'string'   {'mv' 'rv' 'pv' 'pvaf' 'on' 'off'} 'mv';
        'actscale'      'string'   {'on' 'off'}             'off';
        'limcontrib'    'real'     []                       0;
        'topoarg'    'real'     []                       0;
        'gridind'   'real'     []                       0;
        'sumenv'        'string'    {'on' 'off' 'fill'}     'fill'};

    [g varargin] = finputcheck( varargin, fieldlist, 'envtopo_plot', 'ignore');
    if isstr(g), error(g); end;
end

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
%%% Convert g.timerange, g.limits and g.limcontrib to sec from ms %%%%
%
g.timerange = g.timerange/1000;   % the time range of the input data
g.limits(1) = g.limits(1)/1000;   % the time range to plot
if length(g.limits) == 1   % make g.limits at least of length 2
    g.limits(1) = 0; g.limits(2) = 0;
else
    g.limits(2) = g.limits(2)/1000;  %
end;
g.limcontrib = g.limcontrib/1000; % the time range in which to select largest components

%
%%%%%%%%%%%% Collect time range information %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if length(g.limits) == 3 | length(g.limits) > 4 % if g.limits wrong length
    fprintf('envtopo: limits should be 0, [minms maxms], or [minms maxms minuV maxuV].\n');
end

xunitframes = 0; % flag plotting if xmin & xmax are in frames instead of sec
if ~isempty(g.timerange)   % if 'timerange' given
    if g.limits(1)==0 & g.limits(2)==0
        g.limits(1) = min(g.timerange); % if no time 'limits
        g.limits(2) = max(g.timerange); % plot whole data epoch
    end
else % if no 'timerange' given
    if g.limits(1)==0 & g.limits(2)==0 % if no time limits as well,
        fprintf('\nNOTE: No time limits given: using 0 to %d frames\n',frames-1);
        g.limits(1) = 0;
        g.limits(2) = frames-1;
        xunitframes     = 1; % mark frames as time unit instead of sec
    end
end

if isempty(g.timerange)
    xmin = g.limits(1); % (xmin, xmax) are data limits in sec
    xmax = g.limits(2);
else
    xmin = g.timerange(1); % (xmin, xmax) are data limits in sec
    xmax = g.timerange(2);
end

pmin = g.limits(1); % plot min and max sec
if pmin < xmin
    pmin = xmin;     % don't allow plotting beyond the data limits
end
pmax = g.limits(2);
if pmax > xmax
    pmax = xmax;
end

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
%%%%%%%%%%%%%%%%%%%%% Read line color information %%%%%%%%%%%%%%%%%%%%%
%
ENVCOLORS = strvcat('w..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..');

if isempty(g.colorfile)
    g.colorfile = ENVCOLORS; % use default color order above
elseif ~isstr(g.colorfile)
    error('Color file name must be a string.');
end
if strcmpi(g.colorfile,'bold')
    all_bold = 1;
    g.colorfile = ENVCOLORS; % default colors
end
if exist(g.colorfile) == 2  % if an existing file
    cid = fopen(g.colorfile,'r');
    if cid <3,
        error('cannot open color file');
    else
        colors = fscanf(cid,'%s',[3 MAXENVPLOTCHANS]);
        colors = colors';
    end;
else
    colors = g.colorfile;
end
[r c] = size(colors);
for i=1:r
    for j=1:c
        if colors(i,j)=='.',
            if j==1
                error('Color file should have color letter in 1st column.');
            elseif j==2
                colors(i,j)='-';
            elseif j>2
                colors(i,j)=' ';
            end;
        end;
    end;
end;
colors(1,1) = 'k'; % make sure 1st color (for data envelope) is black

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
limitset = 0;
if isempty(g.limits)
    g.limits = 0;
end
if length(g.limits)>1
    limitset = 1;
end

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
    %
    %%%%%% find variance in interval after removing component %%%%%%%%%%%
    %
    if isempty(g.pvaf)
        g.pvaf = g.sortvar;
    end
    if strcmpi(g.pvaf, 'pvaf') | strcmpi('pvaf','on') | strcmpi(g.pvaf,'mv')
        %pvaf(c) = mean(mean((grandERP(:,limframe1:limframe2)-projERP{c}(:,limframe1:limframe2)).^2));
        pvaf(c) = var(reshape(grandERP(:,limframe1:limframe2)-projERP{c}(:,limframe1:limframe2),1,nvals));
    else % if 'pvaf' is 'rv' (or, formerly, 'off')
        %pvaf(c) = mean(mean(projERP{c}(:,limframe1:limframe2).^2));
        pvaf(c) = var(reshape(projERP{c}(:,limframe1:limframe2),1,nvals));
    end;
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
%vardat = mean(mean((grandERP(:,limframe1:limframe2).^2))); % find data variance in interval
vardat = var(reshape(grandERP(:,limframe1:limframe2),1,nvals)); % find data variance in interval

if strcmpi(g.pvaf, 'pvaf') | strcmpi(g.pvaf,'on') | strcmpi(g.pvaf,'mv')
    pvaf = 100-100*pvaf/vardat;
    ot   = g.pvaf;
    if strcmpi(ot,'on'), ot = 'pvaf'; end
elseif strcmpi(g.pvaf, 'rv')| strcmpi(g.pvaf,'off')
    pvaf = 100*pvaf/vardat;
    ot   = 'rv';
end;
[sortpvaf spx] = sort(pvaf);
sortpvaf = sortpvaf(end:-1:1);
spx      = spx(end:-1:1);
npercol = ceil(ncomps/3);

%
%%%%%%%%%%%%%%%%%%%%%%%%% Sort by max variance in data %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(g.pvaf,'mv') | strcmpi(g.pvaf,'on')
    [compvars,compx] = sort(compvars');  % sort clustnums on max variance
else % if 'pvaf' is 'pvaf' or 'rv'
    [tmp,compx]  = sort(pvaf');   % sort clustnums on pvaf or rv
end

compx        = compx(ncomps:-1:1);    % reverse order of sort
compvars     = compvars(ncomps:-1:1)';% reverse order of sort (output var)
compvarorder = g.clustnums(compx);     % actual component numbers (output var)
plotframes   = plotframes(compx);     % plotted comps have these max frames
compframes   = plotframes';           % frame of max variance in each comp (output var)
comptimes    = times(plotframes(compx));  % time of max variance in each comp (output var)
compsplotted = compvarorder(1:ntopos);% (output var)
%
%%%%%%%%%%%%%%%%%%%%%%%% Reduce to ntopos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[plotframes,ifx] = sort(plotframes(1:ntopos));% sort plotframes on their temporal order
plottimes  = times(plotframes);       % convert to times in ms
compx      = compx(ifx);              % indices into clustnums, in plotting order
maporder   = g.clustnums(compx);       % comp. numbers, in plotting order (l->r)
maxproj    = maxproj(:,compx);        % maps in plotting order

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

if strcmpi(g.pvaf,'on') | strcmpi(g.pvaf,'pvaf')
    fprintf('    Component pvaf in interval:  ');
    for t=1:ntopos
        fprintf('%4.2f ',pvaf(t));
    end
    fprintf('\n');
end

sumproj = zeros(size(projERP{1}));
for n = 1:ntopos
    sumproj = sumproj + projERP{compx(n)};
end
varproj = mean(mean((grandERP(:,limframe1:limframe2).^2))); % find data variance in interval
if strcmpi(g.pvaf, 'on')
    sumpvaf = mean(mean((grandERP(:,limframe1:limframe2)-sumproj(:,limframe1:limframe2)).^2));
else
    sumpvaf = mean(mean(sumproj(:,limframe1:limframe2).^2));
end;
if strcmpi(g.pvaf, 'on') | strcmpi(g.pvaf,'pv') | strcmpi(g.pvaf,'pvaf') | strcmpi(g.pvaf,'mv')
    sumpvaf = 100-100*sumpvaf/varproj;
    ot   = 'pvaf';
else % if strcmpi(g.pvaf, 'off') | strcmpi(g.pvaf,'rv')
    sumpvaf = 100*sumpvaf/varproj;
    ot   = 'rv';
end;
if ~xunitframes
    fprintf('    Summed component %s in interval [%4g %4g] ms: %4.2f%%\n',ot, 1000*times(limframe1),1000*times(limframe2), sumpvaf);
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

if isempty(g.limits)
    g.limits = get(axe,'Ylim');
end
set(axe,'GridLineStyle',':')
set(axe,'Xgrid','off')
set(axe,'Ygrid','on')
axes(axe)
set(axe,'Color',axcolor);

%
%%%%%%%%%%%% Collect y-axis range information %%%%%%%%%%%%%%%%%%%%%%%%
%
ylimset = 0; % flag whether hard limits have been set by the user
ymin = min(min(grandERP(:,pframes))); % begin by setting limits from plotted data
ymax = max(max(grandERP(:,pframes)));
if length(g.limits) == 4
    if g.limits(3)~=0 | g.limits(4)~=0 % collect plotting limits from 'limits'
        ymin = g.limits(3);
        ymax = g.limits(4);
        ylimset = 1;
    end
end

%
%%%%%%%%%%%%%%%%% Plot the envelope of the summed selected components %%%%%%%%%%%%%%%%%
%
mapcolors = 1:ntopos+1;

if strcmpi(g.sumenv,'on')  | strcmpi(g.sumenv,'fill')
    sumenv = envelope(sumproj, g.envmode);
    if ~ylimset & max(sumenv) > ymax, ymax = max(curenv); end
    if ~ylimset & min(sumenv) < ymin, ymin = min(curenv); end
    if strcmpi(g.sumenv,'fill')
        %
        % Plot the summed projection filled
        %
        mins = matsel(sumenv,frames,0,2,0);
        p=fill([times times(frames:-1:1)],...
            [matsel(sumenv,frames,0,1,0) mins(frames:-1:1)],FILLCOLOR);
        set(p,'EdgeColor',FILLCOLOR);
        hold on
        %
        % Overplot the data envelope so it is not covered by the fill()'d component
        %
        p=plot(times,matsel(envdata,frames,0,1,1),colors(mapcolors(1),1));% plot the max
        set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)
        p=plot(times,matsel(envdata,frames,0,2,1),colors(mapcolors(1),1));% plot the min
        set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)

    else % if no 'fill'
        tmp = matsel(sumenv,frames,0,2,0);
        p=plot(times,tmp);% plot the min
        hold on
        set(p,'color',FILLCOLOR);
        set(p,'linewidth',2);
        p=plot(times,matsel(sumenv,frames,0,1,0));% plot the max
        set(p,'linewidth',2);
        set(p,'color',FILLCOLOR);
    end
end
if strcmpi(g.pvaf,'pvaf')
    t = text(double(xmin+0.1*(xmax-xmin)), ...
        double(ymin+0.1*(ymax-ymin)), ...
        ['pvaf ' num2str(sumpvaf,'%4.2f') '%']);
    set(t,'fontsize',12,'fontweight','bold')
elseif strcmpi(g.pvaf,'rv') | strcmpi(g.pvaf,'off')
    t = text(double(xmin+0.1*(xmaxxpmin)), ...
        double(ymin+0.1*(ymax-ymin)), ...
        ['rv ' num2str(sumpvaf,'%4.2f') '%']);
    set(t,'fontsize',12,'fontweight','bold')
end

%
% %%%%%%%%%%%%%%%%%%%%%%%% Plot the computed component envelopes %%%%%%%%%%%%%%%%%%
%
envx = [1;compx+1];
for c = 1:ntopos+1
    curenv = matsel(envdata,frames,0,1,envx(c));
    if ~ylimset & max(curenv) > ymax, ymax = max(curenv); end
    p=plot(times,curenv,colors(mapcolors(c),1));% plot the max
    set(gca,'FontSize',12,'FontWeight','Bold')
    if c==1                                % Note: use colors in original
        set(p,'LineWidth',2);              %       component order (if BOLD_COLORS==0)
    else
        set(p,'LineWidth',1);
    end
    if all_bold > 0
        set(p,'LineStyle','-','LineWidth',3);
    elseif mapcolors(c)>15                            % thin/dot 16th-> comp. envs.
        set(p,'LineStyle',':','LineWidth',1);
    elseif mapcolors(c)>10                            %
        set(p,'LineStyle',':','LineWidth',2);
    elseif mapcolors(c)>6                             % dot 6th-> comp. envs.
        set(p,'LineStyle',':','LineWidth',3);
    elseif mapcolors(c)>1
        set(p,'LineStyle',colors(mapcolors(c),2),'LineWidth',1);
        if colors(mapcolors(c),2) == ':'
            set(l1,'LineWidth',2);  % embolden dotted env lines
        end
    end
    hold on
    curenv = matsel(envdata,frames,0,2,envx(c));
    if ~ylimset & min(curenv) < ymin, ymin = min(curenv); end
    p=plot(times,curenv,colors(mapcolors(c),1));% plot the min

    if c==1
        set(p,'LineWidth',2);
    else
        set(p,'LineWidth',1);
    end
    if all_bold > 0
        set(p,'LineStyle','-','LineWidth',3);
    elseif mapcolors(c)>15                            % thin/dot 11th-> comp. envs.
        set(p,'LineStyle',':','LineWidth',1);
    elseif mapcolors(c)>10
        set(p,'LineStyle',':','LineWidth',2);
    elseif mapcolors(c)>6                             % dot 6th-> comp. envs.
        set(p,'LineStyle',':','LineWidth',3);
    elseif mapcolors(c)>1
        set(p,'LineStyle',colors(mapcolors(c),2),'LineWidth',1);
        if colors(mapcolors(c),2) == ':'
            set(l1,'LineWidth',2);  % embolden dotted env lines
        end
    end
    if c==1 & ~isempty(g.vert)
        for v=1:length(g.vert)
            vl=plot([g.vert(v) g.vert(v)], [-1e10 1e10],'k--'); % plot specified vertical lines
            if any(g.limcontrib ~= 0) & v>= length(g.vert)-1;
                set(vl,'linewidth',LIMCONTRIBWEIGHT);
                set(vl,'linestyle',':');
            else
                set(vl,'linewidth',VERTWEIGHT);
                set(vl,'linestyle','--');
            end
        end
    end
    if g.limits(1) <= 0 & g.limits(2) >= 0    % plot vertical line at time zero
        vl=plot([0 0], [-1e10 1e10],'k');
        set(vl,'linewidth',2);
    end

    %
    % plot the n-th component filled
    %
    if g.fillcomp(1)>0 & find(g.fillcomp==c-1)
        fprintf('filling the envelope of component %d\n',c-1);
        mins = matsel(envdata,frames,0,2,envx(c));
        p=fill([times times(frames:-1:1)],...
            [matsel(envdata,frames,0,1,envx(c)) mins(frames:-1:1)],...
            colors(mapcolors(c),1));
        %
        % Overplot the data envlope again so it is not covered by the fill()'d component
        %
        p=plot(times,matsel(envdata,frames,0,1,1),colors(mapcolors(1),1));% plot the max
        set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)
        p=plot(times,matsel(envdata,frames,0,2,1),colors(mapcolors(1),1));% plot the min
        set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)
    end
end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%%%%%%%%%%%%%%%%%%%% Extend y limits by 5% %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~ylimset
    datarange = ymax-ymin;
    ymin = ymin-0.05*datarange;
    ymax = ymax+0.05*datarange;
end
axis([pmin pmax ymin ymax]);

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
    for t=1:ntopos % draw oblique lines from max env vals (or plot top)
        % to map bases, in left to right order
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
        l1 = plot([(plottimes(t)-pmin)/pwidth  ...
            topoleft + 1/pos(3)*(t-1)*1.2*topowidth + (topowidth*0.6)],...
            [data_y 0.68], ...
            colors(linestyles(t)+1)); % 0.68 is bottom of topo maps
        if all_bold > 0
            set(l1,'LineStyle','-','LineWidth',3);
        elseif linestyles(t)>15                        % thin/dot 11th-> comp. envs.
            set(l1,'LineStyle',':','LineWidth',1);
        elseif linestyles(t)>10
            set(l1,'LineStyle',':','LineWidth',2);
        elseif linestyles(t)>5                     % dot 6th-> comp. envs.
            set(l1,'LineStyle',':','LineWidth',3);
        elseif linestyles(t)>1
            set(l1,'LineStyle',colors(linestyles(t)+1,2),'LineWidth',1);
            if colors(linestyles(t)+1,2) == ':'
                set(l1,'LineStyle',colors(linestyles(t)+1,2),'LineWidth',2);
            end
        end
        hold on
        %
        %%%%%%%%%%%%%%%%%%%% add specified vertical lines %%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if g.voffsets(t) > 0
            l2 = plot([(plottimes(t)-xmin)/width  ...
                (plottimes(t)-xmin)/width],...
                [0.6*(maxenv-ymin)/height ...
                0.6*(g.voffsets(t)+maxenv-ymin)/height],...
                colors(linestyles(t)+1));
            if all_bold > 0
                set(l2,'LineStyle','-','LineWidth',3);
            elseif linestyles(t)>15                      % thin/dot 11th-> comp. envs.
                set(l2,'LineStyle',':','LineWidth',1);
            elseif linestyles(t)>10
                set(l2,'LineStyle',':','LineWidth',2);
            elseif linestyles(t)>5                   % dot 6th-> comp. envs.
                set(l2,'LineStyle',':','LineWidth',3);
            else
                set(l1,'LineStyle',colors(linestyles(t)+1,2),'LineWidth',1);
                if colors(linestyles(t)+1,2) == ':'
                    set(l1,'LineWidth',2);
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

    [tmp tmpsort] = sort(maporder);
    [tmp tmpsort] = sort(tmpsort);

    for t=1:ntopos % left to right order  (maporder)
        axt = axes('Units','Normalized','Position',...
            [pos(3)*topoleft+pos(1)+(t-1)*head_sep*topowidth pos(2)+0.66*pos(4) ...
            topowidth topowidth*head_sep]);
        axes(axt)                             % topoplot axes
        cla

        if g.gridind ~= 0
            tmp = zeros(67,67);
            tmp(:)=nan ;
            tmp(g.gridind) = projERP{compx(t)}(:,plotframes(t));
        end
        if ~isempty(varargin)
            figure(myfig);
            toporeplot(tmp, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off',varargin{:});
        else  % if no varargin specified
            figure(myfig);
            toporeplot(tmp, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off','emarkersize',3);
        end

        axis square
        if strcmpi(g.pvaf, 'on') | strcmpi(g.pvaf, 'pvaf') | strcmpi(g.pvaf, 'mv')
            set(gca, 'userdata', ...
                ['text(-0.6, -0.6, ''pvaf: ' sprintf('%6.2f', pvaf(compx(t))) ''');'] );
        else
            set(gca, 'userdata', ...
                ['text(-0.6, -0.6, ''rv: ' sprintf('%6.2f', pvaf(compx(t))) ''');'] );
        end;
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
            if ~isempty(g.clustlabels)
                complabel = g.clustlabels(maporder(t));
            else
                complabel = int2str(maporder(t));        % label comp. numbers
            end
        else
            complabel = compnames(t,:);              % use labels in file
        end
        text(0.00,0.80,complabel,'FontSize',14,...
            'FontWeight','Bold','HorizontalAlignment','Center');
        if strcmpi(g.pvaf, 'on') | strcmpi(g.pvaf, 'pvaf') | strcmpi(g.pvaf, 'mv')
            text(-0.6, -0.6, ['pvaf: ' sprintf('%6.2f', pvaf(compx(t))) ] );
        else
            text(-0.6, -0.6, ['rv: ' sprintf('%6.2f', pvaf(compx(t))) ] );
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

% std_granderp() - compute the grand mean ERP of all the sets in a cluster
%        or of all the sets in the STUDY (based on 'env_erp' input ['contrib'|'all']).
%        This helper function is called by std_envtopo().


function [grandERP, len, STUDY, ALLEEG] = std_granderp(STUDY, ALLEEG, varargin)

if nargin < 4
    help std_granderp;
    return
end

if mod(nargin,2)
    error('std_granderp: Input argument list must be pairs of: ''keyx'', ''valx'' ');
end

clusters = [];
subclus = [];
env_erp = 'contrib';% default grand ERP include only the sets that contribute to the dataset, otherwise 'all' datasets in STUDY.
n = 1;
only_clust = 'off'; % Include only components that were part of the pre-clustering data in the grand ERP

for k = 1:2:(nargin-2)
    switch varargin{k}
        case 'clusters'
            cls = varargin{k+1};
        case 'env_erp'
            env_erp = varargin{k+1};
        case 'only_clust'
            only_clust = varargin{k+1};
        case 'subclus'
            subclus = varargin{k+1};
        case 'baseline'
            baseline = varargin{k+1};
        case 'condition'
            n = varargin{k+1};
    end
end

sets = [];
if strcmpi(env_erp,'all') % all sets in STUDY
    sets = STUDY.setind(n,:);
else % Only sets that contribute to the clusters
    for k = 1:length(cls)
        sets = [sets unique(STUDY.cluster(cls(k)).sets(n,:)) ];
    end
    sets = unique(sets);
end

% Remove clusters from the grand ERP
if exist('subclus')
    if ~isempty(subclus)
        for k = 1: length(sets)
            tmp = [];
            for l = 1:length(subclus)
                for cond = 1:size(STUDY.setind,1)
                    compind = find(STUDY.cluster(subclus(l)).sets(cond,:) == sets(k) );
                    tmp     = [tmp STUDY.cluster(subclus(l)).comps(compind)];
                end;
            end
            subcomps{k} = tmp;
        end
    end
end

fprintf('\n' );
len = length(sets);
for k = 1:len
    EEG = ALLEEG(sets(k));
    if strcmpi(only_clust, 'on')
        comps = STUDY.cluster(cls(1)).preclust.preclustcomps{(mod(sets(k),size(STUDY.setind,2)))}; 
                                                             % Only components that were used in the pre-clustering
    else
        comps = 1:size(EEG.icawinv,2);
    end
    if exist('subcomps')
        if ~isempty(subcomps{k})
            comps = setdiff(comps,subcomps{k});
        end
    end
    fprintf(' .' );
    [tmp_erp] = std_erp(EEG, comps);
    [tmp_scalp] = std_topo(EEG, comps);
    tmp = 0;
    tmp = (tmp_scalp.'*tmp_erp)/len;
    if exist('baseline')
        tmp = rmbase(tmp,EEG.pnts,baseline);
    end
    if k == 1
        grandERP = tmp;
    else
        grandERP = grandERP + tmp;
    end
    ALLEEG(sets(k)) = EEG;
end
