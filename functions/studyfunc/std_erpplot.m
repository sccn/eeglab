% std_erpplot() - Commandline function to plot cluster component ERPs. Either displays 
%                 mean ERP of all requested clusters in the same figure, with ERPs 
%                 for different conditions (if any) plotted in different colors. 
%                 Else, displays ERP for each specified cluster in separate figures 
%                 (per condition), each containing the cluster component ERPs plus 
%                 the grand mean cluster ERP (in bold). ERPs can be plotted only if 
%                 component ERPs were computed and saved in the STUDY EEG datasets. 
%                 These can be computed during pre-clustering using the gui-based 
%                 function pop_preclust() or the equivalent commandline functions 
%                 eeg_createdata() and eeg_preclust(). Called by pop_clustedit().
% Usage:    
%              >> [STUDY] = std_erpplot(STUDY, ALLEEG, key1, val1, key2, val2);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the datasets included 
%                in the STUDY. A STUDY set ALLEEG is typically created by load_ALLEEG().  
%
% Optional inputs:
%   'clusters' - [numeric vector]  -> specific cluster indices to plot.
%                            'all' -> plot all clusters in STUDY {default: 'all'}.
%   'comps'    - [numeric vector]  -> indices of the cluster components to plot.
%                            'all' -> plot all the components in the cluster 
%                {default: 'all'}.
%   'mode'     - ['centroid'|'comps'] plotting mode. 
%                'centroid' -> the average ERPs of the requested clusters are 
%                  plotted in the same figure, with ERPs for  different conditions 
%                  (if any) plotted in different colors. In 'comps' mode, ERPS for 
%                  each specified cluster are plotted in separate figures (per 
%                  condition), each containing cluster component ERPs plus the i
%                  average cluster ERP in bold. Note this parameter has no effect 
%                  if the 'comps' option is used. {default: 'centroid'}.
%   'figure'   - ['on'|'off'] for the 'centroid' mode option. 
%                 'on'  -> plot in a new figure; 
%                 'off' -> plot in the current figure {default: 'on'}
% Outputs:
%   STUDY      - the input STUDY set structure modified with plotted cluster 
%                 mean ERP to allow quick replotting (unless cluster means 
%                 already exists in the STUDY).  
%
%   Example:
%              >> [STUDY] = std_erpplot(STUDY,ALLEEG, 'clusters', 2, 'mode', 'comps');
%                 % Plot cluster-2 components ERPs plus the mean ERP in bold. 
%
%  See also  pop_clustedit(), pop_preclust()
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

% $Log: not supported by cvs2svn $
% Revision 1.10  2006/03/08 21:02:18  arno
% rename func
%
% Revision 1.9  2006/03/08 20:59:23  arno
% rename func
%
% Revision 1.8  2006/03/08 20:59:04  arno
% rename func
%
% Revision 1.7  2006/03/08 20:19:31  arno
% rename func
%
% Revision 1.6  2006/03/07 18:43:13  arno
% plot parent cluster
%
% Revision 1.5  2006/02/16 21:46:29  arno
% header
%
% Revision 1.4  2006/02/16 19:14:49  arno
% now integrating std_plotcompserp.m
%
% Revision 1.3  2006/02/15 22:48:10  arno
% fix rescaling\
%

function STUDY = std_erpplot(STUDY, ALLEEG,  varargin)
icadefs;

% Set default values
cls = []; % plot all clusters in STUDY
mode = 'centroid'; % plot clusters centroid 
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
                    error('std_erpplot: ''clusters'' input takes either specific clusters (numeric vector) or keyword ''all''.');
                end
            end
        case 'comps'
            STUDY = std_plotcomperp(STUDY, ALLEEG,  cls, varargin{k-1});
            return;
        case 'mode' % Plotting mode 'centroid' / 'comps'
            mode = varargin{k-1};
        case 'figure'
            if strcmpi(varargin{k-1},'off')
                figureon = 0;
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
    Ncond = 1;
end
% Plot all the components in the cluster ('comps' mode)
if strcmpi(mode, 'comps')         
    for clus = 1: length(cls) % For each cluster requested
        len = length(STUDY.cluster(cls(clus)).comps);
        if ~isfield(STUDY.cluster(cls(clus)).centroid, 'erp')
            STUDY = std_centroid(STUDY,ALLEEG, cls(clus) , 'erp');
        end
        % For ERP match polarity accross conditions
        for condi = 1: Ncond
            ave_erp(:,condi) = STUDY.cluster(cls(clus)).centroid.erp{condi};
            if  condi == Ncond
                [tmp Avepol] = comppol(ave_erp);
                clear tmp ave_erp
            end
        end
        
        % figure properties
        % -----------------
        figure
        rowcols(2) = ceil(sqrt(Ncond)); 
        rowcols(1) = ceil((Ncond)/rowcols(2));
        pos = get(gcf, 'position');
        magnif = 2.5/sqrt(Ncond);
        set(gcf, 'position', [ pos(1)+15 pos(2)+15 pos(3)*magnif pos(4)/rowcols(2)*rowcols(1)*magnif ]);
        orient tall
        set(gcf,'Color', BACKCOLOR);
        
        for n = 1:Ncond
            try
                clusnval = std_clustread(STUDY, ALLEEG, cls(clus),'erp',n);
            catch,
                warndlg2([ 'Some ERP information is missing, aborting'] , ['Abort - Plot ERP'] );   
                return;
           end
           handl(n) = sbplot(rowcols(1),rowcols(2),n);
           ave_erp = STUDY.cluster(cls(clus)).centroid.erp{n};
           t = STUDY.cluster(cls(clus)).centroid.erp_t;
           [all_erp pol] = comppol(clusnval.erp');
           plot(t/1000,Avepol(n)*all_erp,'color', [0.5 0.5 0.5]);
           hold on
           plot(t/1000,Avepol(n)*ave_erp,'k','linewidth',2);
           xlabel('time [s]');
           ylabel('activations');
           title(['ERP, '  STUDY.cluster(cls(clus)).name ', ' STUDY.condition{n} ', ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) 'Ss']);
           % Make common axis to all conditions
           if n == 1
               ylimits = get(gca,'YLim');
           else
               tmp = get(gca,'YLim');
               ylimits(1) = min(tmp(1),ylimits(1) );
               ylimits(2) = max(tmp(2),ylimits(2) );
           end
           if n == Ncond %set all condition figures to be on the same scale
               for condi = 1: Ncond
                   axes(handl(condi));
                   axis([t(1)/1000 t(end)/1000  ylimits(1)  ylimits(2) ]);
                   axcopy;
               end
           end
        end % finished one condition
    end % finished all requested clusters 
end % Finished 'comps' mode plot option
       
% Plot clusters mean spec/erp
if strcmpi(mode, 'centroid') 
    len = length(cls);
    rowcols(2) = ceil(sqrt(len)); rowcols(1) = ceil((len)/rowcols(2));
    if figureon
        try 
            % optional 'CreateCancelBtn', 'delete(gcbf); error(''USER ABORT'');', 
            h_wait = waitbar(0,['Computing ERP ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,['Computing ERP ...'],'position', [300, 200, 300, 48]);
        end
        figure
    end
    color_codes = {'b', 'r', 'g', 'c', 'm', 'y', 'k','b--', 'r--', 'g--', 'c--', 'm--', 'y--', 'k--','b-.', 'r-.', 'g-.', 'c-.', 'm-.', 'y-.', 'k-.'};
    orient tall
    for k = 1:len % Go through the clusters
        if ~isfield(STUDY.cluster(cls(k)).centroid, 'erp')
            STUDY = std_centroid(STUDY,ALLEEG, cls(k) , 'erp');
        end
        if  (k == 1) 
            erp_min = min(STUDY.cluster(cls(k)).centroid.erp{1});
            erp_max = max(STUDY.cluster(cls(k)).centroid.erp{1});
        end
        if len ~= 1
            handl(k) = sbplot(rowcols(1),rowcols(2),k) ; 
        end
        hold on;
        for n = 1:Ncond
            if k == 1
                leg_color{n} = [STUDY.condition{n}];
            end
            % Compute ERP limits accross conditions and
            % across clusters (all on same scale)  
            erp_min = min(erp_min, min(STUDY.cluster(cls(k)).centroid.erp{n}));
            erp_max = max(erp_max, max(STUDY.cluster(cls(k)).centroid.erp{n}));
            ave_erp(:,n) = STUDY.cluster(cls(k)).centroid.erp{n};
            if n == Ncond
                [ave_erp pol] = comppol(ave_erp);
                t = STUDY.cluster(cls(k)).centroid.erp_t;
                a = [ STUDY.cluster(cls(k)).name ', ' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss'];
                for condi = 1: Ncond
                    plot(t/1000,ave_erp(:,condi),color_codes{condi},'linewidth',2);
                end
            end
            if n == Ncond
                a = [ STUDY.cluster(cls(k)).name ', '  num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ];
                title(a);
                set(gcf,'Color', BACKCOLOR);
                set(gca,'UserData', leg_color);
                set(gcf,'UserData', leg_color);
                if figureon
                    waitbar(k/len,h_wait);
                end
            end
            if (k == len) & (n == Ncond)
                for clsi = 1:len % plot all on same scale
                    if len ~= 1
                        axes(handl(clsi)); 
                    end
                    axis([t(1)/1000 t(end)/1000 erp_min erp_max]);
                    axcopy(gcf, 'leg_color = get(gca,''''UserData'''') ; legend(leg_color); xlabel(''''Time [sec]'''');ylabel(''''Activations'''') ;');
                end
                xlabel('Time [sec]');
                ylabel('Activations');
                if len ~= 1
                    maintitle = ['Average ICA ERP for several clusters across all conditions'];
                    a = textsc(maintitle, 'title'); 
                    set(a, 'fontweight', 'bold'); 
                else
                    a = [ STUDY.cluster(cls(k)).name ' ERP, '  num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ];
                    title(a);
                end
                set(gcf,'Color', BACKCOLOR);
                legend(leg_color);
                if figureon
                    delete(h_wait);
                end
            end
        end % finished the different conditions
    end % finished all clusters 
end % finished 'centroid' plot mode

% std_plotcomperp() - Commandline function, to visualizing cluster component ERPs. 
%                   Displays the ERPs of specified cluster components with the cluster mean 
%                   ERP on separate figures, using one figure for all conditions. 
%                   The ERPs can be visualized only if component ERPs     
%                   were calculated and saved in the EEG datasets in the STUDY.
%                   These can be computed during pre-clustering using the GUI-based function
%                   pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                   and eeg_preclust(). A pop-function that calls this function is pop_clustedit().
% Usage:    
%                   >> [STUDY] = std_plotcomperp(STUDY, ALLEEG, cluster, comps);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   cluster     - single cluster number.  
%
% Optional inputs:
%   comps      - [numeric vector]  -> indices of the cluster components to plot.
%                       'all'                       -> plot all the components in the cluster
%                                                      (as in std_erpplot). {default: 'all'}.
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster
%                     ERP mean, to allow quick replotting (unless cluster mean 
%                     already existed in the STUDY).  
%
%   Example:
%                         >> cluster = 4; comps= 'all';  
%                         >> [STUDY] = std_plotcomperp(STUDY,ALLEEG, cluster, comps);
%                    Plots all components of cluster 4, calls std_erpplot() . 
%
%  See also  pop_clustedit(), pop_preclust()         
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

function STUDY = std_plotcomperp(STUDY, ALLEEG, cls, varargin)
icadefs;

if ~exist('cls')
    error('std_plotcomperp: you must provide a cluster number as an input.');
end
if isempty(cls)
   error('std_plotcomperp: you must provide a cluster number as an input.');
end
if nargin == 3 % no components indices were given
    % Default plot all components of the cluster
    [STUDY] = std_erpplot(STUDY, ALLEEG, 'clusters', cls, 'mode', 'comps');
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
    figure
    orient tall
    set(gcf,'Color', BACKCOLOR);
    for n = 1:Ncond  %for each cond
        abset = STUDY.datasetinfo(STUDY.setind(n,STUDY.cluster(cls).sets(1,comp_ind(ci)))).index;
        subject = STUDY.datasetinfo(STUDY.setind(n,STUDY.cluster(cls).sets(1,comp_ind(ci)))).subject;
        handl(n) = sbplot(rowcols(1),rowcols(2),n);
        hold on
        if Ncond  > 1
            a = [ 'IC' num2str(comp) ' / ' subject ', ' STUDY.cluster(cls).name ', ' STUDY.condition{n} ];
        else
             a = [ 'ERP, IC' num2str(comp) ' / ' subject ', ' STUDY.cluster(cls).name ];
        end
        if ~isfield(STUDY.cluster(cls).centroid, 'erp')
            STUDY = std_centroid(STUDY,ALLEEG, cls, 'erp');
        end
        if ~isfield(ALLEEG(abset).etc,'icaerpparams')
            warndlg2([ 'Dataset ' ALLEEG(abset).filename ' has no ERP info, aborting'] , 'Abort - Plot ERP' ); 
            return;
        end
        [erp, t] = std_readerp(ALLEEG, abset, comp);
        if isempty(erp)
            warndlg2(['eeg_clustedit: file '  ALLEEG(abset).etc.icaerp ' was not found in path ' ALLEEG(abset).filepath], 'Abort - Plot ERP' ); 
            return
        end
        % Change polarity to be similar across conditions
        % (using the ERP centroid).
        if (n == 1) & (Ncond > 1)
            for condi = 1: Ncond
                ave_erp(:,condi) = STUDY.cluster(cls).centroid.erp{condi};
                if  condi == Ncond
                    [tmp Avepol] = comppol(ave_erp);
                    clear tmp ave_erp
                end
            end
        elseif (Ncond == 1)
            Avepol = 1;
        end
        ave_erp = STUDY.cluster(cls).centroid.erp{n};
        [erp tmp] = comppol([erp Avepol(n)*ave_erp]);
        plot(t,erp,'c');
        plot(t,Avepol(n)*ave_erp,'k','linewidth',2);
        xlabel('t [ms]');
        ylabel('activations');
        title(a);
        if n == 1
            ylimits = get(gca,'YLim');
        else
            tmp = get(gca,'YLim');
            ylimits(1) = min(tmp(1),ylimits(1) );
            ylimits(2) = max(tmp(2),ylimits(2) );
        end
        if n == Ncond
            for condi = 1:n
                axes(handl(condi));
                axis([t(1) t(end)  ylimits(1)  ylimits(2) ]);
                axcopy
            end
            if Ncond >1
                textsc([ 'ERP, ' subject ' / IC' num2str(comp) ', ' STUDY.cluster(cls).name ],'title');
            end
        end
    end
end
