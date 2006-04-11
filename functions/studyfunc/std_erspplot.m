% std_erspplot() - plot cluster ERSPs. Displays either mean cluster ERSPs, 
%                  or else all cluster component ERSPs plus the mean cluster 
%                  ERSP in one figure per condition. The ERSPs can be plotted 
%                  only if component ERSPs were computed and saved in the 
%                  EEG datasets in the STUDY. These may either be computed 
%                  during pre-clustering using the gui-based function 
%                  pop_preclust(), or via the equivalent commandline functions 
%                  eeg_createdata() and eeg_preclust(). Called by pop_clustedit().
% Usage:    
%        >> [STUDY] = std_erspplot(STUDY, ALLEEG, key1, val1, key2, val2);  
% Inputs:
%   STUDY    - STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - global vector of EEG structures for the datasets included 
%              in the STUDY. ALLEEG for a STUDY set is typically created 
%              using load_ALLEEG().  
%
% Optional inputs:
%   'clusters' - [numeric vector]  -> cluster numbers to plot.
%                            'all' -> plot all clusters in STUDY 
%                            {default: 'all'}.
%   'comps'    - [numeric vector]  -> cluster components to plot.
%                            'all' -> plot all cluster components 
%                            {default: 'all'}.
%   'channels' - [numeric vector]  -> channels to plot.
%   'mode'     - ['centroid'|'individual'] plotting mode. In 'centroid' 
%                mode, the average ERSPs of the requested clusters or channels  
%                are plotted in the same figure - one per condition. In 
%                'individual' mode, component ERSPs for each
%                cluster (or channel) are plotted in a separate 
%                figure (per condition) with the mean ERSP. 
%                Note that for clusters, this option is irrelevant if component  
%                indices are provided as input {default: 'centroid'} 
%   'figure'  - ['on'|'off'] 'on' -> plot on a new figure; 
%                'off' -> plot on current figure 'figure'.
%                Note: 'off' is optional for one cluster in 'centroid' mode.
%                Useful for incomporating a cluster ERSP into 
%                a complex figure. In case of multiple conditions, 
%                only the first condition is displayed, but clicking on 
%                the figure will open a new figure with all conditions 
%                plotted separately {default: 'on'} 
% Output:
%   STUDY     - the input STUDY set structure modified with plotted cluster 
%               mean ERSPs to allow quick replotting (unless cluster means 
%               already exists in the STUDY).  
%
% Example:
%        >> [STUDY] = std_erspplot(STUDY,ALLEEG, 'clusters', 'all', ...
%                                       'mode', 'centroid');
%           % Plot the mean ERSPs of all clusters in STUDY on the same figure. 
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
% Revision 1.23  2006/03/21 15:43:14  arno
% new .sets format
%
% Revision 1.22  2006/03/20 17:22:45  arno
% new std_clustread
%
% Revision 1.21  2006/03/17 17:54:08  scott
% help msg format
%
% Revision 1.20  2006/03/14 19:57:56  arno
% plottting different conditions
%
% Revision 1.19  2006/03/14 19:56:06  arno
% temporary
% /
%
% Revision 1.18  2006/03/14 19:54:55  arno
% temporary
%
% Revision 1.17  2006/03/12 03:53:36  arno
% fix timevals
%
% Revision 1.16  2006/03/10 15:57:06  arno
% renaming variable
%
% Revision 1.15  2006/03/10 15:53:11  arno
% using log frequencies
%
% Revision 1.14  2006/03/09 23:28:07  arno
% implement new ERSP from Matlab and different structure ec...
%
% Revision 1.13  2006/03/09 22:24:30  arno
% converting frequencies
%
% Revision 1.12  2006/03/09 20:00:24  scott
% help msg
%
% Revision 1.11  2006/03/09 19:05:25  scott
% help msg
%
% Revision 1.10  2006/03/08 20:57:46  arno
% rename func
%
% Revision 1.9  2006/03/08 20:20:12  arno
% rename func
%
% Revision 1.8  2006/03/07 18:43:43  arno
% allow to plot parent cluster
%
% Revision 1.7  2006/02/16 23:42:49  arno
% change tick length
%
% Revision 1.6  2006/02/16 23:33:12  arno
% fix figure off
%
% Revision 1.5  2006/02/16 20:31:39  arno
% integrate std_plotcompersp.m
%
% Revision 1.4  2006/02/16 18:24:03  arno
% title\
%
% Revision 1.3  2006/02/16 00:48:34  arno
% new format etc...
%

function STUDY = std_erspplot(STUDY, ALLEEG,  varargin)
icadefs;

% Set default values
% ------------------
g = finputcheck( varargin, { 'clusters'    { 'integer' 'string' }  { [] [] }           'all';
                             'mode'        'string'     { 'centroid' 'comps' 'individual' } 'centroid';
                             'figure'      'string'     { 'on' 'off' }                 'on';
                             'channels'    'integer'    []                             [];
                             'subject'     'string'     []                             '';
                             'dataset'     'integer'    []                             [];
                             'comps'       'integer'    []                             [] }, 'std_plotersp');
if isstr(g), error(g); end;

% check input parameters
% ----------------------
if isempty(g.clusters)
    g.clusters = 2:length(STUDY.cluster);
end
if strcmpi(g.clusters, 'all')
    g.clusters = 2:length(STUDY.cluster);
end
if ~isfield(STUDY, 'chanmean')
    STUDY.chanmean = [];
end;
Ncond = length(STUDY.condition);

% plot component measures
% -----------------------
if ~isempty(g.comps)
    for ci = 1 : length(g.comps) %for each comp
        if ~isempty(g.clusters)
            abset   = [ STUDY.datasetinfo(STUDY.cluster(g.clusters).sets(:,g.comps(ci))).index ];
            subject =   STUDY.datasetinfo(STUDY.cluster(g.clusters).sets(1,g.comps(ci))).subject;
            comp    =   STUDY.cluster(g.clusters).comps(g.comps(ci));
            fig_title = [ 'ERSP, IC' num2str(comp) ' / ' subject ', ' STUDY.cluster(g.clusters).name ];
            plotcondersp( STUDY, ALLEEG, 'dataset', abset, 'component', comp, 'title', fig_title );
        else
            plotcondersp( STUDY, ALLEEG, 'dataset', g.datasets, 'component', g.comps(ci), 'subject', g.subject);
        end;
    end
    return;
end;

% plot channel measures for one dataset or one subject
% ----------------------------------------------------
if ~isempty(g.channels) & (~isempty(g.dataset) | ~isempty(g.subject))
    for ci = 1 : length(g.channels) %for each comp
        if ~isempty(g.subject), fig_title = [' Subject ' g.subject          ', channel ' int2str(g.channels(ci)) ]; end;
        if ~isempty(g.dataset), fig_title = [' Dataset ' int2str(g.dataset) ', channel ' int2str(g.channels(ci)) ]; end;
        plotcondersp( STUDY, ALLEEG, 'dataset', g.dataset, 'component', g.channels(ci), ...
            'subject', g.subject, 'title', fig_title);
    end
    return;
end;

% -----------------
% compute centroids
% -----------------
if ~isempty(g.clusters) % cluster
    for k = 1:length(g.clusters)
        if ~isfield(STUDY.cluster(g.clusters(k)).centroid,'ersp')
            STUDY = std_centroid(STUDY,ALLEEG, g.clusters(k) , 'ersp');
        end
        if isempty(STUDY.cluster(g.clusters(k)).centroid.ersp)
            warndlg2(['eeg_clustedit: some .icaersp files could not be found for cluster ' STUDY.cluster(g.clusters(k)).name ], 'Abort - Plot ERSP' ); 
            return
        end
    end
else % channel centroid
    for k = 1:length(g.channels)
        STUDY = std_chanmean(STUDY,ALLEEG, g.channels(k), 'ersp');
    end
end;    

% -----------------------------------------------------
% Plot all the components in the cluster ('comps' mode)
% -----------------------------------------------------
if ~isempty(g.clusters) & (strcmpi(g.mode, 'comps') | strcmpi(g.mode, 'individual'))
    for clus = 1: length(g.clusters) % For each cluster requested

        len = length(STUDY.cluster(g.clusters(clus)).comps);
        rowcols(2) = ceil(sqrt(len + 4)); 
        rowcols(1) = ceil((len+4)/rowcols(2));

        % get cluster information
        % -----------------------
        try
            clusncomm = std_clustread(STUDY, ALLEEG, g.clusters(clus),'ersp',[1:Ncond]);
        catch,
            warndlg2([ 'Some ERSP information is missing, aborting'] , ['Abort - Plot ERSP' ] );   
            return;
        end
        
        % plot: one figure per condition
        % ------------------------------
        for n = 1:Ncond
            
            % create new figure
            % -----------------
            figure
            orient tall
            maintitle = ['ERSP, cond. ' num2str(n) ', ' STUDY.cluster(g.clusters(clus)).name ];
            a = textsc(maintitle,'title'); 
            set(a, 'fontweight', 'bold'); 
            set(gcf,'Color', BACKCOLOR);
            sbplot(rowcols(1),rowcols(2),[1 rowcols(2)+2 ]);
            
            % plot centroid
            % -------------
            ave_ersp = STUDY.cluster(g.clusters(clus)).centroid.ersp{n};
            idat     = STUDY.datasetinfo(STUDY.cluster(g.clusters(clus)).sets(1,1)).index;
            ersp_times = STUDY.cluster(g.clusters(clus)).centroid.ersp_times;
            ersp_freqs = STUDY.cluster(g.clusters(clus)).centroid.ersp_freqs;
            fig_title  = [ STUDY.cluster(g.clusters(clus)).name ' average ERSP, ' ...
                num2str(length(unique(STUDY.cluster(g.clusters(clus)).sets(1,:)))) 'Ss' ];
            %lim = STUDY.cluster(g.clusters(clus)).centroid.ersp_limits{n}; %plotting limits
            plotersp( ave_ersp, ersp_times, ersp_freqs, 'title', fig_title);

            % Plot the individual component ersp 
            % ----------------------------------
            for k = 1:len %  
                abset   = STUDY.datasetinfo(STUDY.cluster(g.clusters(clus)).sets(n,k)).index;
                subject = STUDY.datasetinfo(STUDY.cluster(g.clusters(clus)).sets(n,k)).subject;
                comp    = STUDY.cluster(g.clusters(clus)).comps(k);
                                
                fig_title = [ 'ic' num2str(comp) '/' subject ];
                if k <= rowcols(2) - 2 %first sbplot row
                    sbplot(rowcols(1),rowcols(2),k+2); 
                else  %other sbplot rows
                    sbplot(rowcols(1),rowcols(2),k+4);  
                end
                
                plotersp( clusncomm.ersp{n,k}(:,:), clusncomm.ersp_times, clusncomm.ersp_freqs, ...
                    'title', fig_title, 'xlabel', 'off', 'ylabel', 'off', 'cbar', 'off');
                set(get(gca,'Title'),'FontSize',8);
                drawnow;
            end % finished all components in cluster 
       end % finished all conditions
    end % finished all requested clusters 
end % Finished 'comps' mode plot option

% -----------------------------------------
% Plot all the channels in the mean channel
% -----------------------------------------
if ~isempty(g.channels) & strcmpi(g.mode, 'individual')
    for chan = 1: length(g.channels) % For each cluster requested

        len = size(STUDY.setind,2);
        rowcols(2) = ceil(sqrt(len + 4)); 
        rowcols(1) = ceil((len+4)/rowcols(2));

        % get cluster information
        % -----------------------
        try
            channcomn = std_chanread(STUDY, ALLEEG, g.channels(chan),'ersp',[1:Ncond]);
        catch,
            warndlg2([ 'Some ERSP information is missing, aborting'] , ['Abort - Plot ERSP' ] );   
            return;
        end
        
        % plot: one figure per condition
        % ------------------------------
        for n = 1:Ncond
            
            % create new figure
            % -----------------
            figure
            orient tall
            maintitle = ['ERSP, cond. ' num2str(n) ', Channel ' int2str(g.channels) ];
            a = textsc(maintitle,'title'); 
            set(a, 'fontweight', 'bold'); 
            set(gcf,'Color', BACKCOLOR);
            sbplot(rowcols(1),rowcols(2),[1 rowcols(2)+2 ]);
            
            % plot centroid
            % -------------
            ave_ersp   = STUDY.chanmean.ersp{n, g.channels(chan)};
            ersp_times = STUDY.chanmean.ersp_times;
            ersp_freqs = STUDY.chanmean.ersp_freqs;
            fig_title  = [ 'Channel ' int2str(g.channels(chan)) ' average ERSP, ' num2str(length(STUDY.datasetinfo)) 'Ss' ];
            %lim = STUDY.cluster(g.clusters(clus)).centroid.ersp_limits{n}; %plotting limits
            plotersp( ave_ersp, ersp_times, ersp_freqs, 'title', fig_title);

            % Plot the individual component ersp 
            % ----------------------------------
            for k = 1:len %  
                abset   = STUDY.datasetinfo(STUDY.setind(n,k)).index;
                subject = STUDY.datasetinfo(STUDY.setind(n,k)).subject;
                                
                fig_title = [ 'chan ' num2str(g.channels(chan)) '/' subject ];
                if k <= rowcols(2) - 2 %first sbplot row
                    sbplot(rowcols(1),rowcols(2),k+2); 
                else  %other sbplot rows
                    sbplot(rowcols(1),rowcols(2),k+4);  
                end
                
                plotersp( channcomn.ersp{n,k}, channcomn.ersp_times, channcomn.ersp_freqs, ...
                    'title', fig_title, 'xlabel', 'off', 'ylabel', 'off', 'cbar', 'off');
                set(get(gca,'Title'),'FontSize',8);
                drawnow;
            end % finished all components in cluster 
       end % finished all conditions
    end % finished all requested clusters 
end % Finished 'comps' mode plot option

% ------------------------------
% Plot only the cluster centroid
% ------------------------------
if strcmpi(g.mode, 'centroid') & ~isempty(g.clusters) 
    
    % create figure
    % -------------
    len = length(g.clusters);
    if len ~= 1
        rowcols(2) = ceil(sqrt(len)/Ncond); rowcols(1) = ceil((len)/rowcols(2)); rowcols(2) = rowcols(2)*Ncond;
    else
        rowcols(2) = Ncond;  rowcols(1) = 1;
    end
    if strcmpi(g.figure, 'on')
        figure
        pos = get(gcf, 'position');
        magnif = 2.5/sqrt(Ncond);
        set(gcf, 'position', [ pos(1)+15 pos(2)+15 pos(3)*magnif pos(4)/rowcols(2)*rowcols(1)*magnif ]);
        orient tall
        set(gcf,'Color', BACKCOLOR);
    end;    
    if length(g.clusters) > 1
        maintitle = ['Average ERSP for clusters ' int2str(g.clusters(:)') ];
        a = textsc(maintitle, 'title'); 
        set(a, 'fontweight', 'bold'); 
    end;
    
    for k = 1:len 
        
        % find maximum value
        % ------------------
        maxval = 0;
        for n = 1:Ncond
            maxval = max(max(max(abs(STUDY.cluster(g.clusters(k)).centroid.ersp{n}))), maxval);
        end;
 
        if strcmpi(g.figure, 'off') % only for summary mode: average all conditions
            % average all conditions
            % ----------------------
            for n = 1:Ncond          
                if n == 1
                    ave_ersp = STUDY.cluster(g.clusters(k)).centroid.ersp{n}/Ncond;
                else
                    ave_ersp = ave_ersp + STUDY.cluster(g.clusters(k)).centroid.ersp{n}/Ncond;
                end;
            end;
            logfreqs = STUDY.cluster(g.clusters(k)).centroid.ersp_freqs;
            timevals = STUDY.cluster(g.clusters(k)).centroid.ersp_times;
            plotersp( ave_ersp,timevals,logfreqs, 'clim', [-maxval maxval], 'title', 'Average ERSP');
        else
            
            for n = 1:Ncond
                sbplot(rowcols(1),rowcols(2),(k-1)*Ncond+n), 
                a = [ STUDY.cluster(g.clusters(k)).name ' ERSP, ' num2str(length(unique(STUDY.cluster(g.clusters(k)).sets(1,:)))) 'Ss, ' STUDY.condition{n}];
                ave_ersp = STUDY.cluster(g.clusters(k)).centroid.ersp{n};
                logfreqs = STUDY.cluster(g.clusters(k)).centroid.ersp_freqs;
                timevals = STUDY.cluster(g.clusters(k)).centroid.ersp_times;
                plotersp( ave_ersp,timevals,logfreqs, 'clim', [-maxval maxval], 'title', a, 'cbar', 'off');                
                if (k-1)*Ncond+n > (rowcols(1)-1)*rowcols(2)
                    xlabel('Time [ms]');
                else
                    xlabel('');
                end;
                if n ~= 1
                    set(gca,'ytick',[]);
                    set(gca,'yticklabel', []);
                    ylabel('');
                end;
                if n == Ncond
                    cbar;
                end;
            end;
        end % Finish plotting all centroids for one condition
    end  % Finished all conditions
end % Finished 'centroid' mode plot option

% ------------------------------
% Plot only the cluster centroid
% ------------------------------
if strcmpi(g.mode, 'centroid') & ~isempty(g.channels) 
    
    % create figure
    % -------------
    len = length(g.channels);
    if len ~= 1
        rowcols(2) = ceil(sqrt(len)/Ncond); rowcols(1) = ceil((len)/rowcols(2)); rowcols(2) = rowcols(2)*Ncond;
    else
        rowcols(2) = Ncond;  rowcols(1) = 1;
    end
    if strcmpi(g.figure, 'on')
        figure
        pos = get(gcf, 'position');
        magnif = 2.5/sqrt(Ncond);
        set(gcf, 'position', [ pos(1)+15 pos(2)+15 pos(3)*magnif pos(4)/rowcols(2)*rowcols(1)*magnif ]);
        orient tall
        set(gcf,'Color', BACKCOLOR);
    end;    
    if length(g.channels) > 1
        maintitle = ['Average ERSP for channels ' int2str(g.channels(:)') ];
        a = textsc(maintitle, 'title'); 
        set(a, 'fontweight', 'bold'); 
    end;
    
    for k = 1:len 
        
        % find maximum value
        % ------------------
        maxval = 0;
        for n = 1:Ncond
            maxval = max(max(max(abs(STUDY.chanmean.ersp{n, g.channels(k)}))), maxval);
        end;
 
        if strcmpi(g.figure, 'off') % only for summary mode: average all conditions
            % average all conditions
            % ----------------------
            for n = 1:Ncond          
                if n == 1
                    ave_ersp = STUDY.chanmean.ersp{n, g.channels(k)}/Ncond;
                else
                    ave_ersp = ave_ersp + STUDY.chanmean.ersp{n, g.channels(k)}/Ncond;
                end;
            end;
            logfreqs = STUDY.chanmean.ersp_freqs;
            timevals = STUDY.chanmean.ersp_times;
            plotersp( ave_ersp,timevals,logfreqs, 'clim', [-maxval maxval], 'title', 'Average ERSP');
        else
            
            for n = 1:Ncond
                sbplot(rowcols(1),rowcols(2),(k-1)*Ncond+n), 
                a = [ 'Channel ' int2str(g.channels(k)) ' ERSP, ' num2str(length(STUDY.subject)) 'Ss, ' STUDY.condition{n}];
                ave_ersp = STUDY.chanmean.ersp{n, g.channels(k)};
                logfreqs = STUDY.chanmean.ersp_freqs;
                timevals = STUDY.chanmean.ersp_times;
                plotersp( ave_ersp,timevals,logfreqs, 'clim', [-maxval maxval], 'title', a, 'cbar', 'off');                
                if (k-1)*Ncond+n > (rowcols(1)-1)*rowcols(2)
                    xlabel('Time [ms]');
                else
                    xlabel('');
                end;
                if n ~= 1
                    set(gca,'ytick',[]);
                    set(gca,'yticklabel', []);
                    ylabel('');
                end;
                if n == Ncond
                    cbar;
                end;
            end;
        end % Finish plotting all centroids for one condition
    end  % Finished all conditions
end % Finished 'centroid' mode plot option

% ------------------------------------------------------------
% std_plotcompersp() - plot ERSP for one or several conditions
% ------------------------------------------------------------
function STUDY = plotcondersp(STUDY, ALLEEG, varargin)

    g = finputcheck( varargin, { 'dataset'     'integer'     []                  [];
                                 'component'   'integer'     []                  [];
                                 'condition'   'string'      []                  '';
                                 'subject'     'string'      []                  '';
                                 'channel'     'integer'     []                  [];
                                 'figure'      'string'      { 'on' 'off' }      'on';
                                 'title'       'string'      []                  '';
                                 'rowcols'     'integer'     []                  [1 1];
                                 'rowcolpos'   'integer'     []                  1;
                                 'plotargs'    'cell'        []                  {}}, 'plotcondersp');
    if isstr(g), error(g); end;

    % conditions
    % ----------
    Ncond = length(STUDY.condition);
    if isempty(g.dataset)
        if isempty(g.condition), g.condition = [1:Ncond];
        else                     g.condition = strmatch(g.condition, lower(STUDY.condition));
        end;
    else
        g.condition = 1;
    end;
    if isempty(g.condition), error('Unknown condition'); end;
    Ncond = length(g.condition);
        
    % figure properties
    % -----------------
    if strcmpi(g.figure, 'on')
        icadefs;
        figure
        pos = get(gcf, 'position');
        magnif = 2.5/sqrt(Ncond); if Ncond == 1, magnif = 1; end;
        g.rowcols(2) = ceil(sqrt(Ncond)); 
        g.rowcols(1) = ceil((Ncond)/g.rowcols(2));
        set(gcf, 'position', [ pos(1)+15 pos(2)+15 pos(3)*magnif pos(4)/g.rowcols(2)*g.rowcols(1)*magnif ]);
        orient tall
        set(gcf,'Color', BACKCOLOR);
        g.rowcolpos = 1;
    end;
   
    % dataset indices
    % ---------------
    if ~isempty(g.subject)
       g.dataset = STUDY.setind(:,strmatch(lower(g.subject), lower(STUDY.subject)));
       if isempty(g.dataset), error('Could not find subject'); end;
    end;
    subject = STUDY.datasetinfo(g.dataset(1)).subject;

    % retrieve data
    % -------------
    if ~isempty(g.component)
        [ersp, logfreqs, timeval] = std_readersp(ALLEEG, g.dataset, g.component, STUDY.preclust.erspclusttimes,  STUDY.preclust.erspclustfreqs);
    else
        [ersp, logfreqs, timeval] = std_readersp(ALLEEG, g.dataset, -g.channel, STUDY.preclust.erspclusttimes,  STUDY.preclust.erspclustfreqs);
    end;    

    % find index in setind
    % --------------------
    for n = g.condition  %for each cond
        
        sbplot(g.rowcols(1),g.rowcols(2),g.rowcolpos);
        g.rowcolpos = g.rowcolpos + 1;
        idat = g.dataset(n);
        if Ncond >1, fig_title = [ g.title ', ' STUDY.condition{n} ]; 
        else         fig_title =   g.title; end;
        
        plotersp( ersp(:,:,n), timeval, logfreqs, 'clim', [-4 4], 'title', fig_title, g.plotargs{:});
    end

% sub function to plot ERSPs
% --------------------------
function plotersp( data, times, freqs, varargin);

    g = finputcheck(varargin, { 'title'    'string'    []               '';
                                'xlabel'   'string'    { 'on' 'off' }   'on';
                                'ylabel'   'string'    { 'on' 'off' }   'on';
                                'cbar'     'string'    { 'on' 'off' }   'on';
                                'clim'     'real'      []               [] });
    if isstr(g), error(g); end;
    logfreqs = log(freqs);

    tftopo( data, times,logfreqs,'limits', [times(1) times(end) logfreqs(1) logfreqs(end) g.clim ],...
        'title', g.title, 'verbose', 'off', 'axcopy', 'off');
    ft = str2num(get(gca,'yticklabel'));
    ft = exp(1).^ft;
    ft = unique(round(ft));
    ftick = get(gca,'ytick');
    ftick = exp(1).^ftick;
    ftick = unique(round(ftick));
    ftick = log(ftick);
    set(gca,'ytick',ftick);
    set(gca,'yticklabel', num2str(ft));
    if strcmpi(g.xlabel, 'on')
        xlabel('Time [ms]');
    else    
        xlabel('');
        set(gca, 'xtick', []);
    end;
    if strcmpi(g.ylabel, 'on')
        ylabel('Power (dB)');
    else    
        ylabel('');
        set(gca, 'ytick', []);
    end;
    
    axcopy(gcf, [' ft = str2num(get(gca,''''yticklabel'''')); ft = exp(1).^ft; ft = unique(round(ft)); fti = get(gca,''''ytick''''); fti = exp(1).^fti; fti = unique(round(fti));'...
        'fti = log(fti); set(gca, ''''ytick'''',fti); set(gca, ''''yticklabel'''',num2str(ft)); xlabel(''''Time [ms]''''); ylabel(''''Frequency [Hz]''''); cbar; clear ft fti;' ]);
    if strcmpi(g.cbar, 'on')
        cbar;
    end;