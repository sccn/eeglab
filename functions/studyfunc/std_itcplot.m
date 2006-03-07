% cls_plotclustitc() - Commandline function, to visualizing cluster/s ITCs. 
%                   Displays either mean cluster/s ITC/s, or all cluster/s component
%                   ITCs with the mean cluster/s ITC in one figure (per cluster & condition).
%                   The ITCs can be visualized only if component ITCs     
%                   were calculated and saved in the EEG datasets in the STUDY.
%                   These can be computed during pre-clustering using the GUI-based function
%                   pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                   and eeg_preclust(). A pop-function that calls this function is pop_clustedit().
% Usage:    
%                   >> [STUDY] = cls_plotclustitc(STUDY, ALLEEG, key1, val1, key2, val2);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Optional inputs:
%   'clusters'   - [numeric vector]  -> specific cluster numbers to plot.
%                     'all'                         -> plot all clusters in STUDY.
%                     {default: 'all'}.
%   'comps'      - [numeric vector]  -> indices of the cluster components to plot.
%                       'all'        -> plot all the components in the cluster {default: 'all'}.
%   'mode'       - ['centroid'|'comps'] a plotting mode. In 'centroid' mode, the average ERSPs 
%                     of the requested clusters are plotted in the same figure - one per condition. 
%                     In 'comps' mode, component ERSPs for each cluster are plotted in a
%                     separate figure (per condition) with the cluster mean ERSP.
%                     {default: 'centroid'}. Note that this option is irrelevant if component
%                     indices are provided as input.
%   'figure'   - ['on'|'off'] plots on a new figure ('on')  or plots on current
%                    figure ('off'). 'figure' 'off' is optional for one cluster in 'centroid' mode.
%                    Useful for incomporating cluster ITC into a complex figure.
%                    In case of multiple conditions only the first condition is displayed,
%                    once plotted clicking on the figure will open a new figure with 
%                    all the conditions plotted. {default: 'on'}. 
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster 
%                     mean ITCs, to allow quick replotting (unless cluster means 
%                     already exists in the STUDY).  
%
%   Example:
%                         >> [STUDY] = cls_plotclustitc(STUDY,ALLEEG, 'clusters', 'all', 'mode', 'centroid');
%                    Plots the mean ITCs of all the clusters in STUDY on the same figure. 
%
%  See also  pop_clustedit, pop_preclust, eeg_createdata, cls_plotcompitc         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

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

function STUDY = cls_plotclustitc(STUDY, ALLEEG,  varargin)
icadefs;
% Set default values
cls = []; % plot all clusters in STUDY
mode = 'centroid'; % plot clusters centroid 
figureon = 1;

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
                    error('cls_plotclustitc: ''clusters'' input takes either specific clusters (numeric vector) or keyword ''all''.');
                end
            end
        case 'comps'
            STUDY = cls_plotcompitc(STUDY, ALLEEG,  cls, varargin{k-1});
            return;
        case 'mode' % Plotting mode 'centroid' / 'comps'
            mode = varargin{k-1};
        case 'figure' % plot on exisiting figure or a new figure
            if strcmpi(varargin{k-1},'off')
                if length(cls) == 1 % only in case of one cluster
                    figureon = 0;
                end
            end
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

Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond =1;
end
% Plot all the components in the cluster ('comps' mode)
if strcmpi(mode, 'comps')         
    for clus = 1: length(cls) % For each cluster requested
        len = length(STUDY.cluster(cls(clus)).comps);
        try 
            h_wait = waitbar(0,['Plotting ITC ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,['Plotting ITC ...'],'position', [300, 200, 300, 48]);
        end
        rowcols(2) = ceil(sqrt(len + 4)); rowcols(1) = ceil((len+4)/rowcols(2));
        if ~isfield(STUDY.cluster(cls(clus)).centroid,'itc')
            STUDY = cls_centroid(STUDY,ALLEEG, cls(clus) , 'itc');
        end
        for n = 1:Ncond
            try
                clusncomm = cls_clusread(STUDY, ALLEEG, cls(clus),'itc',n);
            catch,
                warndlg2([ 'Some ITC information is missing, aborting'] , ['Abort - Plot ITC' ] );   
                delete(h_wait)
                return;
           end
           figure
           orient tall
            maintitle = ['ITC, cond. ' num2str(n) ', ' STUDY.cluster(cls(clus)).name ];
            a = textsc(maintitle,'title'); 
            set(a, 'fontweight', 'bold'); 
            set(gcf,'Color', BACKCOLOR);
            sbplot(rowcols(1),rowcols(2),[1 rowcols(2)+2 ]) ,
            ave_itc = STUDY.cluster(cls(clus)).centroid.itc{n};
            %lim = STUDY.cluster(cls(clus)).centroid.itc_limits{n}; %plotting limits
            itc_times = ALLEEG(STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(cls(clus)).sets(1,1))).index).etc.icaerspparams.times;
            logfreqs = STUDY.cluster(cls(clus)).centroid.itc_logf;
            a = [ STUDY.cluster(cls(clus)).name ' average ITC, ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) 'Ss' ];
            tftopo(abs(ave_itc),itc_times,logfreqs,'limits', [itc_times(1) itc_times(end) logfreqs(1) logfreqs(end) -.5 .5],...
                'title', a, 'verbose', 'off', 'axcopy', 'off');
            ft = str2num(get(gca,'yticklabel'));
            ft = exp(1).^ft;
            ft = unique(round(ft));
            ftick = get(gca,'ytick');
            ftick = exp(1).^ftick;
            ftick = unique(round(ftick));
            ftick = log(ftick);
            set(gca,'ytick',ftick);
            set(gca,'yticklabel', num2str(ft));
            xlabel('Time [ms]');
            axcopy(gcf, [' ft = str2num(get(gca,''''yticklabel'''')); ft = exp(1).^ft; ft = unique(round(ft)); fti = get(gca,''''ytick''''); fti = exp(1).^fti; fti = unique(round(fti));'...
                'fti = log(fti); set(gca, ''''ytick'''',fti); set(gca, ''''yticklabel'''',num2str(ft)); xlabel(''''Time [ms]''''); cbar(''pos''); clear ft fti;' ]);
            cbar('pos');
            for k = 1:len % Plot the individual component itc  
                abset = STUDY.datasetinfo(STUDY.setind(n,STUDY.cluster(cls(clus)).sets(1,k))).index;
                subject = STUDY.datasetinfo(STUDY.setind(n,STUDY.cluster(cls(clus)).sets(1,k))).subject;
                comp = STUDY.cluster(cls(clus)).comps(k);
                params = ALLEEG(abset).etc.icaitcparams;
                a = [ 'ic' num2str(comp) '/' subject];
                if k <= rowcols(2) - 2 %first sbplot row
                    sbplot(rowcols(1),rowcols(2),k+2); 
                else  %other sbplot rows
                    sbplot(rowcols(1),rowcols(2),k+4);  
                end
                tftopo(abs(clusncomm.itc{k}),params.times,clusncomm.logf{k},'limits', ...
                    [params.times(1) params.times(end) clusncomm.logf{k}(1) clusncomm.logf{k}(end) -.5 .5],...
                    'title', a, 'verbose', 'off', 'axcopy', 'off');
                set(gca, 'xtick', [], 'ytick', []);
                set(get(gca,'Title'),'FontSize',8)
                xlabel('');
                ylabel('');
                axcopy(gcf, [' ft = str2num(get(gca,''''yticklabel'''')); ft = exp(1).^ft; ft = unique(round(ft)); fti = get(gca,''''ytick''''); fti = exp(1).^fti; fti = unique(round(fti));'...
                'fti = log(fti); set(gca, ''''ytick'''',fti); set(gca, ''''yticklabel'''',num2str(ft)); xlabel(''''Time [ms]''''); cbar(''pos''); clear ft fti;' ]);
                waitbar((k*n)/(Ncond*len),h_wait);
                if k == len
                    cbar('pos');
                end
            end % finished all components in cluster 
       end % finished all conditions
       delete(h_wait);
    end % finished all requested clusters 
end % Finished 'comps' mode plot option

% Plot clusters mean itc
if strcmpi(mode, 'centroid') 
    len = length(cls);
    if len ~= 1
        rowcols(2) = ceil(sqrt(len)/Ncond); rowcols(1) = ceil((len)/rowcols(2)); rowcols(2) = rowcols(2)*Ncond;
    else
        rowcols(2) = Ncond;  rowcols(1) = 1;
    end
    if figureon
        try 
            % optional 'CreateCancelBtn', 'delete(gcbf); error(''USER ABORT'');', 
            h_wait = waitbar(0,['Computing ITC ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,['Computing ITC ...'],'position', [300, 200, 300, 48]);
        end
    end
    
    % Compute cluster centroid
    % ------------------------
    for k = 1:len 
        if ~isfield(STUDY.cluster(cls(k)).centroid,'itc')
            STUDY = cls_centroid(STUDY,ALLEEG, cls(k) , 'itc');
        end
        if isempty(STUDY.cluster(cls(k)).centroid.itc)
            warndlg2(['eeg_clustedit: some .icaitc files could not be found for cluster ' STUDY.cluster(cls(k)).name ], 'Abort - Plot itc' ); 
            return
        end
        % ITC plotting limitis is the average limits across all clusters
        %for n = 1:Ncond
        %    lim(n) = lim(n)+STUDY.cluster(cls(k)).centroid.itc_limits{n}; %plotting limits
        %end
    end

    params = ALLEEG(STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(cls(1)).sets(1,1))).index).etc.icaerspparams;
    if figureon
        figure
        pos = get(gcf, 'position');
        magnif = 2.5/sqrt(Ncond);
        set(gcf, 'position', [ pos(1)+15 pos(2)+15 pos(3)*magnif pos(4)/rowcols(2)*rowcols(1)*magnif ]);
        orient tall
        set(gcf,'Color', BACKCOLOR);
    end;
    
    if len > 1
        maintitle = ['Average ERSP for clusters ' int2str(cls(1:len)) ];
        a = textsc(maintitle, 'title'); 
        set(a, 'fontweight', 'bold'); 
    end;
    
    for k = 1:len 
        
        % find maximum value
        % ------------------
        maxval = 0;
        for n = 1:Ncond
            maxval = max(max(max(abs(STUDY.cluster(cls(k)).centroid.itc{n}))), maxval);
        end;
        
        % plot
        % ----
       if ~figureon % only for summary mode: average all conditions
            % average all conditions
            % ----------------------
            for n = 1:Ncond          
                if n == 1
                    ave_itc = STUDY.cluster(cls(k)).centroid.itc{n}/Ncond;
                else
                    ave_itc = ave_itc + STUDY.cluster(cls(k)).centroid.itc{n}/Ncond;
                end;
            end;
            logfreqs = STUDY.cluster(cls(k)).centroid.itc_logf;
            tftopo(abs(ave_itc),params.times,logfreqs,'limits', [params.times(1) params.times(end) logfreqs(1) logfreqs(end) -maxval maxval],...
                   'title', 'Average ITC', 'verbose', 'off');
            ft = str2num(get(gca,'yticklabel'));
            ft = exp(1).^ft;
            ft = unique(round(ft));
            ftick = get(gca,'ytick');
            ftick = exp(1).^ftick;
            ftick = unique(round(ftick));
            ftick = log(ftick);
            set(gca,'ytick',ftick);
            set(gca,'yticklabel', num2str(ft));
            xlabel('Time [ms]');
            cbar('pos');
       else
           for n = 1:Ncond
               sbplot(rowcols(1),rowcols(2),(k-1)*Ncond+n), 
               a = [ STUDY.cluster(cls(k)).name ' ITC, ' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss, ' STUDY.condition{n}];
               ave_itc  = STUDY.cluster(cls(k)).centroid.itc{n};
               logfreqs = STUDY.cluster(cls(k)).centroid.itc_logf;
               tftopo(abs(ave_itc),params.times,logfreqs,'limits', [params.times(1) params.times(end) logfreqs(1) logfreqs(end) -maxval maxval],...
                      'title', a, 'verbose', 'off');
               ft = str2num(get(gca,'yticklabel'));
               ft = exp(1).^ft;
               ft = unique(round(ft));
               ftick = get(gca,'ytick');
               ftick = exp(1).^ftick;
               ftick = unique(round(ftick));
               ftick = log(ftick);
               set(gca,'ytick',ftick);
               set(gca,'yticklabel', num2str(ft));
               if (k-1)*Ncond+n > (rowcols(1)-1)*rowcols(2)
                   xlabel('Time [ms]');
               else
                   xlabel('');
               end;
               cbar('pos');
               waitbar((k*n)/(len*Ncond),h_wait);
           end;
        end % Finish plotting all centroids for one condition
    end  % Finished all conditions
    if figureon
        delete(h_wait)
    end
end % Finished 'centroid' mode plot option

% cls_plotcompitc() - Commandline function, to visualizing cluster component ITC images. 
%                    Displays the ITC images of specified cluster components on separate figures,
%                    using one figure for all conditions. 
%                   The ITCs can be visualized only if component ITCs     
%                   were calculated and saved in the EEG datasets in the STUDY.
%                   These can be computed during pre-clustering using the GUI-based function
%                   pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                   and eeg_preclust(). A pop-function that calls this function is pop_clustedit().
% Usage:    
%                   >> [STUDY] = cls_plotcompitc(STUDY, ALLEEG, cluster, comps);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   cluster     - single cluster number.  
%
% Optional inputs:
%   comps      - [numeric vector]  -> indices of the cluster components to plot.
%                       'all'                       -> plot all the components in the cluster
%                                                      (as in cls_plotclustmap). {default: 'all'}.
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster itc
%                     image mean, to allow quick replotting (unless cluster mean 
%                     already existed in the STUDY).  
%
%   Example:
%                         >> cluster = 4; comps= [1 7 10];  
%                         >> [STUDY] = cls_plotcompmap(STUDY,ALLEEG, cluster, comps);
%                    Plots components 1, 7 & 10  itcs of cluster 4 on separate figures. 
%
%  See also  pop_clustedit, pop_preclust, eeg_createdata, cls_plotclustitc         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $ Log: cls_plotcompitc.m,v $

function STUDY = cls_plotcompitc(STUDY, ALLEEG, cls, varargin)
icadefs;
if ~exist('cls')
    error('cls_plotcompitc: you must provide a cluster number as an input.');
end
if isempty(cls)
   error('cls_plotcompitc: you must provide a cluster number as an input.');
end
if nargin == 3 % no components indices were given
    % Default plot all components of the cluster
    [STUDY] = cls_plotclustitc(STUDY, ALLEEG, 'clusters', cls, 'mode', 'comps');
    return
else
    comp_ind = varargin{1}; 
end

Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond =1;
end
for ci = 1 : length(comp_ind) %for each comp
   rowcols(2) = ceil(sqrt(Ncond)); rowcols(1) = ceil((Ncond)/rowcols(2));
   comp = STUDY.cluster(cls).comps(comp_ind(ci));     
   subject = STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(cls).sets(1,comp_ind(ci)))).subject;
   
   % figure properties
   % -----------------
   figure
   pos = get(gcf, 'position');
   magnif = 2.5/sqrt(Ncond);
   set(gcf, 'position', [ pos(1)+15 pos(2)+15 pos(3)*magnif pos(4)/rowcols(2)*rowcols(1)*magnif ]);
   orient tall
   set(gcf,'Color', BACKCOLOR);
   
   for n = 1:Ncond  %for each cond
        abset = STUDY.datasetinfo(STUDY.setind(n,STUDY.cluster(cls).sets(1,comp_ind(ci)))).index;
        if ~isfield(ALLEEG(abset).etc,'icaitcparams')
            warndlg2([ 'Dataset ' num2str(abset) ' has no ITC info, aborting'] , ['Abort - Plot ITC']); 
            return;
        end
        params = ALLEEG(abset).etc.icaitcparams;
        sbplot(rowcols(1),rowcols(2),n), 
        [itc, logfreqs] = cls_readitc(ALLEEG, abset, comp);
        if isempty(itc)
            warndlg2(['eeg_clustedit: file '  ALLEEG(abset).etc.icalogitc ' was not found in path ' ALLEEG(abset).filepath], 'Abort - Plot ITC' ); 
            return
        end
        if Ncond >1
            a = [ 'ITC, IC' num2str(comp) ' / ' subject ', ' STUDY.cluster(cls).name ', ' STUDY.condition{n} ];
        else
            a = ['ITC, IC' num2str(comp) ' / ' subject  ', ' STUDY.cluster(cls(clus)).name];
        end
        tftopo(abs(itc),params.times,logfreqs,'limits', [params.times(1) params.times(end) logfreqs(1) logfreqs(end) -.5 .5],...
            'title', a, 'verbose', 'off', 'axcopy', 'off');
                ft = str2num(get(gca,'yticklabel'));
        ft = exp(1).^ft;
        ft = unique(round(ft));
        ftick = get(gca,'ytick');
        ftick = exp(1).^ftick;
        ftick = unique(round(ftick));
        ftick = log(ftick);
        set(gca,'ytick',ftick);
        set(gca,'yticklabel', num2str(ft));
        cbar('pos');
        axcopy(gcf, [' ft = str2num(get(gca,''''yticklabel'''')); ft = exp(1).^ft; ft = unique(round(ft)); fti = get(gca,''''ytick''''); fti = exp(1).^fti; fti = unique(round(fti));'...
            'fti = log(fti); set(gca, ''''ytick'''',fti); set(gca, ''''yticklabel'''',num2str(ft)); xlabel(''''Time [ms]''''); cbar(''pos''); clear ft fti;' ]);
   end
end
