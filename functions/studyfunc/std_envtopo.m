% cls_envtopo() - Creates envtopo image for a STUDY set, uses clusters contribution instead of components.  
%             Plot the envelope of a data epoch, plus envelopes and average scalp maps of specified 
%             or largest-contributing clusters for each condition. Click on individual axes to examine them in detail. 
% Usage:
%             >> cls_envtopo(STUDY, ALLEEG, 'key1', 'val1', ...);
% Inputs:
%   ALLEEG       - ALLEEG data structure, can also be one EEG dataset structure.
%   STUDY        - an eeglab STUDY set (that contains some of these EEGs structures)
%
% Optional inputs:
%  'clusters'  = [integer array] vector of cluster numbers. If one cluster plots the 
%                  cluster contribution to the data envelope. If multiple
%                  clusters selects the largest contributing clusters from within 
%                  and plots their contribution. {default| [] -> 'all'}
%  'subclus'  = [integer array] vector of cluster numbers to remove from data before  
%                  computing the ERP envelope of the data.   
%  'env_erp'   = ['contrib' | 'all'] 
%                  'contrib' - > for one cluster contribution to the grand ERP envelope,
%                       the envelope includes only the datasets that are part of that cluster.
%                  'all' -> grand ERP envelope includes all datasets in STUDY, in case of  
%                  multiple clusters that is the only possible option. 
%  'only_clust' = [ 'on' | 'off'] Which dataset components to include in the grand ERP. 
%                     'on' - include only the components that were part of the clustering 
%                        (for example, if components fromt the clustering because of high rv 
%                        don't include their contribution in the grand ERP).
%                     'off' - include all the components of the datasets 
%                        (execpt the ones in the subtructed clusters). {default 'off'}. 
%  'baseline'    = [minms maxms] - a new baseline to remove from the grand
%                       ERP and cluster ERP contributions.
%  'diff'           = [condition1 condition2] the indexes of two condition.
%                        Plots additional figure with the difference of the
%                        two condition envtopo. 
%  'clustnums'  = [integer array] vector of cluster numbers to plot {default|0 -> all}
%                  Else if int < 0, the number of largest contributing clusters to plot 
%                  {default|[] -> 7}
%  'timerange' = start and end input data latencies (in ms) {default: from 'limits' if any}
%  'limits'    = 0 or [minms maxms] or [minms maxms minuV maxuV]. Specify start/end plot
%                  (x) limits (in ms) and min/max y-axis limits (in uV). If 0, or if both
%                  minmx & maxms == 0 -> use latencies from 'timerange' (else 0:frames-1).
%                  If both minuV and maxuV == 0 -> use data uV limits {default: 0}
%  'limcontrib' = [minms maxms]  time range (in ms) in which to rank clusters contribution
%                  (boundaries shown with thin dotted lines) {default|[]|[0 0] -> plotting limits}
%  'vert'       = vector of times (in ms) at which to plot vertical dashed lines {default|[] -> none}
%
% see also: envtopo, cls_granderp, envtopo_plot
%

function cls_envtopo(STUDY, ALLEEG, varargin)
icadefs;
if nargin < 4
    help STUDY_envtopo;
    return
end

if mod(nargin,2)
    error('STUDY_envtopo: Input argument list must be pairs of: ''keyx'', ''valx'' ');
end

clusters = [];
subclus = []; 
e_options{1} = 'env_erp';
e_options{2} = 'contrib';% default grand ERP include only the sets that contribute to the dataset, otherwise 'all' datasets in STUDY.
e_options{3} = 'only_clust'; % Include only components that were part of the pre-clustering data in the grand ERP
e_options{4} = 'off'; % Default is off, so include all components of the sets (execpt ones in the subtructed cluster) 
p_options = {};

for k = 1:2:(nargin-2)
    switch varargin{k}
        case 'clusters'
            clusters = varargin{k+1};
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
            t = ALLEEG(STUDY.datasetinfo(STUDY.setind(1)).index).times; % convert baseline from ms to indexes
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

if isempty(clusters) | strcmpi(clusters,'all') % Default all clusters in STUDY
    clusters = [];
    for k = 2:length(STUDY.cluster)
        if ~strncmpi('Notclust',STUDY.cluster(k).name,8)
            clusters = [clusters k];
        end
    end
end
% If some of the requestes clusters to subtruct are in clusters, remove them  
clusters = setdiff(clusters,subclus);
%if length(clusters) > 1
%    e_options{2} = 'all';
%end
Ncond = length(STUDY.condition); %number of conditions
if Ncond == 0
    Ncond = 1;
end 

for n = 1:Ncond
    % compute the grand mean ERP envelope of the cluster (for specific
   % condition).
   fprintf('\n Computing grand ERP for condition %d.', n); 
   [grandERP, set_len, STUDY, ALLEEG] = cls_granderp(STUDY, ALLEEG, 'clusters', clusters, 'condition', n, e_options{:});
   grandERPtot{n} = grandERP;
   if ~exist('limits') | (length(limits)  == 2)%set limits to be same for all conditions
       tmpmin = min(min(grandERP)); 
       tmpmax = max(max(grandERP)); 
       datarange = tmpmax-tmpmin;
       tmpmin = tmpmin-0.05*datarange;
       tmpmax = tmpmax+0.05*datarange;
       if n == 1 
           ymin = tmpmin; 
           ymax = tmpmax; 
       else
           ymin = min(tmpmin,ymin); 
           ymax = max(tmpmax,ymax); 
       end
       if n == Ncond
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

for n = 1:Ncond
    % compute the grand mean ERP envelope of the cluster (for specific
   % condition).
	for cls = 1:length(clusters)
       len = length(STUDY.cluster(clusters(cls)).comps);
        try
            clusscalp = cls_clusread(STUDY, ALLEEG, clusters(cls),'scalp');
        catch,
            warndlg2([ 'Some topoplot information is missing, aborting'] , 'Abort - STUDY_envtopo' );   
            return;
       end
       try
           cluserp = cls_clusread(STUDY, ALLEEG, clusters(cls),'erp', n);
           if exist('baseline')
               cluserp.erp = rmbase(cluserp.erp,ALLEEG(STUDY.datasetinfo(STUDY.setind(1)).index).pnts,baseline);
           end
       catch,
            warndlg2([ 'Some ERP information is missing, aborting'] , 'Abort - STUDY_envtopo' );   
            return;
       end   
       % compute grand mean back projection ERP for the cluster
       projERP = 0;
       fprintf('\n Computing projected component ERP of cluster %d: ', (clusters(cls)) ); 
       val_ind = find(~isnan(clusscalp.grid{1}(:))); % find non-NAN values
       for k = 1:len
           tmp = clusscalp.grid{k}(val_ind);
           projERP = projERP + tmp*cluserp.erp(k,:);
           fprintf('.'); 
       end
       tot_projERP{cls} = projERP/set_len;
       clus_names{cls} = [STUDY.cluster(clusters(cls)).name];
       numind = strfind(STUDY.cluster(clusters(cls)).name,' ')+1; %find the cluster number
       clus_ind{cls} = [STUDY.cluster(clusters(cls)).name(numind(end):end) ', '];
   end
   p_options{end+1} = 'clustlabels';
   p_options{end+1} = clus_names;
   p_options{end+1} = 'gridind';
   p_options{end+1} = val_ind;
   p_options{end+1} = 'sortvar';
   p_options{end+1} = 'pvaf';
   % compute the grand mean ERP envelope of the cluster (for specific
   % condition).   
   if cls == 1
       figure; set(gcf,'Color', BACKCOLOR); 
       orient landscape;
       % for memory saving once a computation of a condition is done the
       % grand ERP for it is cleared grandERPtot(1), which make the current
       % condition grand ERP always in index 1
       envtopo_plot(grandERPtot{1},tot_projERP,'envmode' ,'avg', 'fillcomp', 1, ...
           'dispmaps', 'on', 'title', ['Projected cluster: ' clus_ind{:} ' ' STUDY.condition{n}], p_options{:} );
   else
       figure; set(gcf,'Color', BACKCOLOR);
       orient landscape;
       envtopo_plot(grandERPtot{1},tot_projERP,'envmode' ,'avg', ...
           'dispmaps', 'on', 'title', [ STUDY.condition{n}], p_options{:} );
   end
   %chil = get(gcf, 'Children');
   %childf(n) = chil(end);
    %if n == 1
    %    ylimits = get(childf(n),'YLim');
    %else
    %    tmp = get(childf(n),'YLim');
    %    ylimits(1) = min(tmp(1),ylimits(1) );
    %    ylimits(2) = max(tmp(2),ylimits(2) );
    %end
   if exist('diffc')
      if  ~exist('diff_h')
          diff_h = figure;
          orient landscape;
          set(gcf,'Color', BACKCOLOR);
          if length(clus_ind) < 5
              maintitle =  ['Projected clusters: ' clus_ind{:}];
          else
              maintitle =  ['Projected clusters difference between conditions'];
          end
          a = textsc(maintitle, 'title'); 
          set(a, 'fontweight', 'bold'); 
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
               diff_proj{cls}  = proj1{1} - proj2{1};
               proj1(1) = []; %saving memory
               proj2(1) =[];
           end
           clear proj1 proj2
           %plot difference
           subplot(3,1,3)
%           if exist('limits')
%               p_options{find(strcmpi({p_options{:}},'limits'))+1} = [limits(1:2) limits(3:4)/2];
%           end
           envtopo_plot(erp1-erp2,diff_proj,'envmode' ,'avg', 'dispmaps','off','title','Difference',p_options{:} );
%           tmp =  get(gcf, 'Children');
%           set(tmp(4),'YLim',ylimits);
%           set(tmp(6),'YLim',ylimits);
           figure
           if cls == 1
               envtopo_plot(erp1-erp2,diff_proj,'envmode' ,'avg', 'fillcomp', 1, ...
                    'dispmaps', 'on', 'title', 'Difference between the two conditions', p_options{:} );
            else
                envtopo_plot(erp1-erp2,diff_proj,'envmode' ,'avg', 'dispmaps','on','title','Difference between the two conditions',p_options{:} );
            end
           orient landscape;
       end
   end
   grandERPtot(1) = [];
end
