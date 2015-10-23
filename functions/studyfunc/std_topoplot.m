% std_topoplot() - Command line function to plot cluster component and mean scalp maps. 
%                  Displays either mean cluster/s scalp map/s, or all cluster/s components
%                  scalp maps with the mean cluster/s scsalp map in one figure.
%                  The scalp maps can be visualized only if component scalp maps     
%                  were calculated and saved in the EEG datasets in the STUDY.
%                  These can be computed during pre-clustering using the GUI-based function
%                  pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                  and eeg_preclust(). A pop-function that calls this function is 
%                  pop_clustedit().
% Usage:    
%              >> [STUDY] = std_topoplot(STUDY, ALLEEG, key1, val1, key2, val2);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in 
%                the STUDY. ALLEEG for a STUDY set is typically created using load_ALLEEG().
% Optional inputs:
%   'clusters'   - [numeric vector| 'all']  -> specific cluster numbers to plot.
%                  'all' -> plot all clusters in STUDY.
%                  {default: 'all'}.
%   'comps'      - [numeric vector | 'all']  -> indices of the cluster components to plot.
%                  'all' -> plot all the components in the cluster 
%                  {default: 'all'}.
%   'mode'       - ['together'|'apart'] a plotting mode. In 'together' mode, the average 
%                  scalp maps of the requested clusters are plotted in the same figure,
%                  one per condition.  In 'apart' mode, component scalp maps for each 
%                  cluster are plotted in a separate figure for each condition, plus the 
%                  cluster mean map. Note that this option is irrelevant if component 
%                  indices ('comps' above) are provided.{default: 'apart'}. 
%   'figure'     - ['on'|'off'] for the 'together' mode option, plots on
%                  a new figure ('on'),  or on the current figure ('off').
%                  {default: 'on'}.
% Outputs:
%   STUDY        - the input STUDY set structure modified with plotted cluster scalp
%                  map means, to allow quick replotting (unless clusters meands 
%                  already exists in th STUDY).  
%
%   Example:
%              % Plot the mean scalp maps for clusters 1 through 20 on the same figure. 
%           >> [STUDY] = std_topoplot(STUDY,ALLEEG, 'clusters', [1:20], 'mode', 'together');
%
%  See also  pop_clustedit(), pop_preclust()
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 07, 2005, hilit@sccn.ucsd.edu
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

function STUDY = std_topoplot(STUDY, ALLEEG,  varargin)
icadefs;

% Set default values
cls = [];
mode = 'together'; % plot all cluster mean scalp maps on one figure
figureon = 1; % plot on a new figure

for k = 3:2:nargin
    switch varargin{k-2}
        case 'clusters'
            if isnumeric(varargin{k-1})
                cls = varargin{k-1};
                if isempty(cls)
                    cls = 2:length(STUDY.cluster);
                end
            else
                if isstr(varargin{k-1}) & strcmpi(varargin{k-1}, 'all')
                    cls = 2:length(STUDY.cluster);
                else
                    error('std_topoplot: ''clusters'' input takes either specific clusters (numeric vector) or keyword ''all''.');
                end
            end
        case 'plotsubjects' % legacy
            mode = 'apart';
        case 'comps'
            if isstr( varargin{k-1} ), mode = 'apart'; 
            else STUDY = std_plotcompmap(STUDY, ALLEEG,  cls, varargin{k-1}); return;
            end;
        case 'mode' % Plotting mode 'together' / 'apart'
            mode = varargin{k-1};
        case 'figure'
            if strcmpi(varargin{k-1},'off')
                figureon = 0;
            end
        case 'plotrad'
            inputPlotrad = varargin{k-1};
    end
end
        
% select clusters to plot
% -----------------------
if isempty(cls)
    tmp =[];
    cls = 2:length(STUDY.cluster); % plot all clusters in STUDY
    for k = 1: length(cls)
        % don't include 'Notclust' clusters
        if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) & ~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)
            tmp = [tmp cls(k)];
        end
    end
    cls = tmp;
end;

% Plot all the components in the cluster
disp('Drawing components of cluster (all at once)...');
if ~isfield(STUDY.cluster,'topo'), STUDY.cluster(1).topo = []; end;
for clus = 1: length(cls) % For each cluster requested
    if isempty(STUDY.cluster(cls(clus)).topo)
        STUDY = std_readtopoclust(STUDY,ALLEEG, cls(clus));
    end;
end
if strcmpi(mode, 'apart')         
    for clus = 1: length(cls) % For each cluster requested
        len = length(STUDY.cluster(cls(clus)).comps);
        if len > 0 % A non-empty cluster 
            h_topo = figure;
            rowcols(2) = ceil(sqrt(len + 4)); rowcols(1) = ceil((len+4)/rowcols(2));
            clusscalp = STUDY.cluster(cls(clus));
            ave_grid = clusscalp.topo;
            tmp_ave = ave_grid;
            tmp_ave(find(isnan(tmp_ave))) = 0; % remove NaN values from grid for later correlation calculation.  
            for k = 1:len
                abset   = STUDY.datasetinfo(STUDY.cluster(cls(clus)).sets(1,k)).index;
                subject = STUDY.datasetinfo(STUDY.cluster(cls(clus)).sets(1,k)).subject;
                comp = STUDY.cluster(cls(clus)).comps(k);
                [Xi,Yi] = meshgrid(clusscalp.topoy,clusscalp.topox);                     
                scalpmap = squeeze(clusscalp.topoall{k}); % already correct polarity
                if k <= rowcols(2) - 2 %first sbplot row
                    figure(h_topo);
                    sbplot(rowcols(1),rowcols(2),k+2) , 
                    toporeplot(scalpmap, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off','xsurface', Xi, 'ysurface', Yi );
                    title([subject '/' 'ic' num2str(comp)   ], 'interpreter', 'none');
                    colormap(DEFAULT_COLORMAP);
                    %waitbar(k/(len+1),h_wait)
                else %other sbplot rows
                    figure(h_topo)
                    sbplot(rowcols(1),rowcols(2),k+4) , 
                    toporeplot(scalpmap, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off','xsurface', Xi, 'ysurface', Yi );
                    title([subject '/' 'ic' num2str(comp)], 'interpreter', 'none');
                    colormap(DEFAULT_COLORMAP);
                    %waitbar(k/(len+1),h_wait)
                end
            end
            figure(h_topo)
            sbplot(rowcols(1),rowcols(2),[1 rowcols(2)+2 ]) ,
            toporeplot(ave_grid, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off');
            title([ STUDY.cluster(cls(clus)).name ' (' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) ' Ss, ' num2str(length(STUDY.cluster(cls(clus)).comps)),' ICs)']);
            %title([ STUDY.cluster(cls(clus)).name ' average scalp map, ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) 'Ss']);
            set(gcf,'Color', BACKCOLOR);
            colormap(DEFAULT_COLORMAP);
            %waitbar(1,h_wait)
            %delete(h_wait)
            orient tall  % fill the figure page for printing
            axcopy
        end % Finished one cluster plot 
    end   % Finished plotting all clusters
end % Finished 'apart' plotting mode

% Plot clusters centroid maps
if strcmpi(mode, 'together') 
    len = length(cls);
    rowcols(2) = ceil(sqrt(len)); rowcols(1) = ceil((len)/rowcols(2));
    if figureon
        try 
            % optional 'CreateCancelBtn', 'delete(gcbf); error(''USER ABORT'');', 
            h_wait = waitbar(0,'Computing topoplot ...', 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,'Computing topoplot ...','position', [300, 200, 300, 48]);
        end
        figure   
    end
    for k = 1:len 
       if len ~= 1
            sbplot(rowcols(1),rowcols(2),k)  
       end
        tmpcmap = colormap(DEFAULT_COLORMAP);
        toporeplot(STUDY.cluster(cls(k)).topo, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off','colormap', tmpcmap);
        title([ STUDY.cluster(cls(k)).name ' (' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) ' Ss, ' num2str(length(STUDY.cluster(cls(k)).comps)),' ICs)']);
        colormap(DEFAULT_COLORMAP);
        %title([ STUDY.cluster(cls(k)).name ', ' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ]);
        if figureon
            waitbar(k/len,h_wait)
        end
    end
    if figureon
        delete(h_wait)
    end
    if len ~= 1
        maintitle = 'Average scalp map for all clusters';
        a = textsc(maintitle, 'title'); 
        set(a, 'fontweight', 'bold');
        set(gcf,'name', maintitle);
    else
        title([ STUDY.cluster(cls(k)).name ' (' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) ' Ss, ' num2str(length(STUDY.cluster(cls(k)).comps)),' ICs)']);
        set(gcf,'name',['Scalp map of ' STUDY.cluster(cls(k)).name ' (' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) ' Ss, ' num2str(length(STUDY.cluster(cls(k)).comps)),' ICs)']);
        %title([ STUDY.cluster(cls(k)).name ' scalp map, ' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ]);
    end
    set(gcf,'Color', BACKCOLOR);

    orient tall
    axcopy
end        

% std_plotcompmap() - Commandline function, to visualizing cluster components scalp maps. 
%                   Displays the scalp maps of specified cluster components on separate figures. 
%                   The scalp maps can be visualized only if component scalp maps     
%                   were calculated and saved in the EEG datasets in the STUDY.
%                   These can be computed during pre-clustering using the GUI-based function
%                   pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                   and eeg_preclust(). A pop-function that calls this function is pop_clustedit().
% Usage:    
%                   >> [STUDY] = std_plotcompmap(STUDY, ALLEEG, cluster, comps);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   cluster     - single cluster number.  
%
% Optional inputs:
%   comps      - [numeric vector]  -> indices of the cluster components to plot.
%                       'all'                       -> plot all the components in the cluster
%                                                      (as in std_topoplot). {default: 'all'}.
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster scalp
%                     map mean, to allow quick replotting (unless cluster mean 
%                     already existed in the STUDY).  
%
%   Example:
%                         >> cluster = 4; comps= [1 7 10];  
%                         >> [STUDY] = std_plotcompmap(STUDY,ALLEEG, cluster, comps);
%                    Plots components 1, 7 & 10  scalp maps of cluster 4 on separate figures. 
%
%  See also  pop_clustedit, pop_preclust, eeg_createdata, std_topoplot         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 07, 2005, hilit@sccn.ucsd.edu
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

function STUDY = std_plotcompmap(STUDY, ALLEEG, cls, varargin)
icadefs;

if ~exist('cls')
    error('std_plotcompmap: you must provide a cluster numberas an input.');
end
if isempty(cls)
   error('std_plotcompmap: you must provide a cluster numberas an input.');
end
if nargin == 3 % no component indices were given
    % Default: plot all components of the cluster
    [STUDY] = std_topoplot(STUDY, ALLEEG, 'clusters', cls, 'mode', 'apart');
    return
else
    comp_ind = varargin{1}; 
end
STUDY = std_readtopoclust(STUDY,ALLEEG, cls);
for ci = 1:length(comp_ind)
    abset = STUDY.datasetinfo(STUDY.cluster(cls).sets(1,comp_ind(ci))).index;
    subject = STUDY.datasetinfo(STUDY.cluster(cls).sets(1,comp_ind(ci))).subject;
    comp  = STUDY.cluster(cls).comps(comp_ind(ci));
    grid  = STUDY.cluster(cls).topoall{comp_ind(ci)};
    xi    = STUDY.cluster(cls).topox;
    yi    = STUDY.cluster(cls).topoy;
    [Xi,Yi] = meshgrid(yi,xi);
    figure;
    toporeplot(grid, 'style', 'both', 'plotrad',0.5,'intrad',0.5,'xsurface', Xi, 'ysurface', Yi, 'verbose', 'off');
    title([subject ' / ' 'IC' num2str(comp) ', ' STUDY.cluster(cls).name ]);
    set(gcf,'Color', BACKCOLOR);
    colormap(DEFAULT_COLORMAP);
    axcopy;
end

