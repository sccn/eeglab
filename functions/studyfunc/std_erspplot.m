% cls_plotclustersp() - Commandline function, to visualizing cluster/s ERSPs. 
%                   Displays either mean cluster/s ERSP/s, or all cluster/s component 
%                   ERSPs with the mean cluster/s ERSP in one figure (per condition).
%                   The ERSPs can be visualized only if component ERSPs     
%                   were calculated and saved in the EEG datasets in the STUDY.
%                   These can be computed during pre-clustering using the GUI-based function
%                   pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                   and eeg_preclust(). A pop-function that calls this function is pop_clustedit().
% Usage:    
%                   >> [STUDY] = cls_plotclustersp(STUDY, ALLEEG, key1, val1, key2, val2);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Optional inputs:
%   'clusters'   - [numeric vector]  -> specific cluster numbers to plot.
%                     'all'          -> plot all clusters in STUDY.
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
%                    Useful for incomporating cluster ERSP into a complex figure.
%                    In case of multiple conditions only the first condition is displayed,
%                    once plotted clicking on the figure will open a new figure with 
%                    all the conditions plotted. {default: 'on'}. 
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster 
%                     mean ERSPs, to allow quick replotting (unless cluster means 
%                     already exists in the STUDY).  
%
%   Example:
%                         >> [STUDY] = cls_plotclustersp(STUDY,ALLEEG, 'clusters', 'all', 'mode', 'centroid');
%                    Plots the mean ERSPs of all the clusters in STUDY on the same figure. 
%
%  See also  pop_clustedit, pop_preclust, eeg_createdata, cls_plotcompersp         
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
% Revision 1.6  2006/02/16 23:33:12  arno
% fix figure off
%
% Revision 1.5  2006/02/16 20:31:39  arno
% integrate cls_plotcompersp.m
%
% Revision 1.4  2006/02/16 18:24:03  arno
% title\
%
% Revision 1.3  2006/02/16 00:48:34  arno
% new format etc...
%

function STUDY = cls_plotclustersp(STUDY, ALLEEG,  varargin)
icadefs;

% Set default values
cls = 2:length(STUDY.cluster); % plot all clusters in STUDY
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
                    error('cls_plotclustersp: ''clusters'' input takes either specific clusters (numeric vector) or keyword ''all''.');
                end
            end
        case 'comps'
            STUDY = cls_plotcompersp(STUDY, ALLEEG,  cls, varargin{k-1});
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

tmp =[];
for k = 1: length(cls)
    % don't include 'Notclust' clusters
    if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) & ~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)
        tmp = [tmp cls(k)];
    end
end
cls = tmp;
clear tmp

Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond = 1;
end
% Plot all the components in the cluster ('comps' mode)
if strcmpi(mode, 'comps')         
    for clus = 1: length(cls) % For each cluster requested
        len = length(STUDY.cluster(cls(clus)).comps);
        try 
            h_wait = waitbar(0,['Plotting ERSP ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,['Plotting ERSP ...'],'position', [300, 200, 300, 48]);
        end
        rowcols(2) = ceil(sqrt(len + 4)); rowcols(1) = ceil((len+4)/rowcols(2));
        if ~isfield(STUDY.cluster(cls(clus)).centroid,'ersp')
            STUDY = cls_centroid(STUDY,ALLEEG, cls(clus) , 'ersp');
        end
        try
            clusncomm = cls_clusread(STUDY, ALLEEG, cls(clus),'ersp',[1:Ncond]);
        catch,
            warndlg2([ 'Some ERSP information is missing, aborting'] , ['Abort - Plot ERSP' ] );   
            delete(h_wait)
            return;
       end
       for n = 1:Ncond
           figure
           orient tall
            maintitle = ['ERSP, cond. ' num2str(n) ', ' STUDY.cluster(cls(clus)).name ];
            a = textsc(maintitle,'title'); 
            set(a, 'fontweight', 'bold'); 
            set(gcf,'Color', BACKCOLOR);
            sbplot(rowcols(1),rowcols(2),[1 rowcols(2)+2 ]) ,
            ave_ersp = STUDY.cluster(cls(clus)).centroid.ersp{n};
            lim = STUDY.cluster(cls(clus)).centroid.ersp_limits{n}; %plotting limits
            ersp_times = ALLEEG(STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(cls(clus)).sets(1,1))).index).etc.icaerspparams.times;
            logfreqs = STUDY.cluster(cls(clus)).centroid.ersp_logf;
            a = [ STUDY.cluster(cls(clus)).name ' average ERSP, ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) 'Ss' ];
            tftopo(ave_ersp,ersp_times,logfreqs,'limits', [ersp_times(1) ersp_times(end) logfreqs(1) logfreqs(end)],...
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
                'fti = log(fti); set(gca, ''''ytick'''',fti); set(gca, ''''yticklabel'''',num2str(ft)); xlabel(''''Time [ms]''''); ylabel(''''Frequency [Hz]''''); cbar; clear ft fti;' ]);
            cbar;
            for k = 1:len % Plot the individual component ersp  
                abset = STUDY.datasetinfo(STUDY.setind(n,STUDY.cluster(cls(clus)).sets(1,k))).index;
                subject = STUDY.datasetinfo(STUDY.setind(n,STUDY.cluster(cls(clus)).sets(1,k))).subject;
                comp = STUDY.cluster(cls(clus)).comps(k);
                params = ALLEEG(abset).etc.icaerspparams;
                a = [ 'ic' num2str(comp) '/' subject ];
                if k <= rowcols(2) - 2 %first sbplot row
                    sbplot(rowcols(1),rowcols(2),k+2); 
                else  %other sbplot rows
                    sbplot(rowcols(1),rowcols(2),k+4);  
                end
                tftopo(clusncomm.ersp{k}(:,:,n),params.times,clusncomm.logf{k},'limits', ...
                    [params.times(1) params.times(end) clusncomm.logf{k}(1) clusncomm.logf{k}(end) -lim lim],...
                    'title', a, 'verbose', 'off', 'axcopy', 'off');
                set(gca, 'xtick', [], 'ytick', []);
                set(get(gca,'Title'),'FontSize',8)
                xlabel('');
                ylabel('');
                axcopy(gcf, [' ft = str2num(get(gca,''''yticklabel'''')); ft = exp(1).^ft; ft = unique(round(ft)); fti = get(gca,''''ytick''''); fti = exp(1).^fti; fti = unique(round(fti));'...
                     'fti = log(fti); set(gca, ''''ytick'''',fti); set(gca, ''''yticklabel'''',num2str(ft)); xlabel(''''Time [ms]''''); ylabel(''''Frequency [Hz]''''); cbar; clear ft fti;' ]);
                waitbar((k*n)/(Ncond*len),h_wait);
                if k == len
                    cbar;
                end
            end % finished all components in cluster 
       end % finished all conditions
       delete(h_wait);
    end % finished all requested clusters 
end % Finished 'comps' mode plot option

% Plot clusters mean ersp
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
            h_wait = waitbar(0,['Computing ERSP ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,['Computing ERSP ...'],'position', [300, 200, 300, 48]);
        end
    end
    % ERSP plotting limitis is the average limits across all clusters

    % Compute cluster centroid
    % ------------------------
    for k = 1:len 
        if ~isfield(STUDY.cluster(cls(k)).centroid,'ersp')
            STUDY = cls_centroid(STUDY,ALLEEG, cls(k) , 'ersp');
        end
        if isempty(STUDY.cluster(cls(k)).centroid.ersp)
            warndlg2(['eeg_clustedit: some .icaersp files could not be found for cluster ' STUDY.cluster(cls(k)).name ], 'Abort - Plot ERSP' ); 
            return
        end
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
            maxval = max(max(max(abs(STUDY.cluster(cls(k)).centroid.ersp{n}))), maxval);
        end;
 
        if ~figureon % only for summary mode: average all conditions
            % average all conditions
            % ----------------------
            for n = 1:Ncond          
                if n == 1
                    ave_ersp = STUDY.cluster(cls(k)).centroid.ersp{n}/Ncond;
                else
                    ave_ersp = ave_ersp + STUDY.cluster(cls(k)).centroid.ersp{n}/Ncond;
                end;
            end;
            logfreqs = STUDY.cluster(cls(k)).centroid.ersp_logf;
            tftopo(ave_ersp,params.times,logfreqs,'limits', [params.times(1) params.times(end) logfreqs(1) logfreqs(end) -maxval maxval],...
                   'title', 'Average ERSP', 'verbose', 'off');
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
            cbar;
        else
            for n = 1:Ncond
                sbplot(rowcols(1),rowcols(2),(k-1)*Ncond+n), 
                a = [ STUDY.cluster(cls(k)).name ' ERSP, ' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss, ' STUDY.condition{n}];
                ave_ersp = STUDY.cluster(cls(k)).centroid.ersp{n};
                logfreqs = STUDY.cluster(cls(k)).centroid.ersp_logf;
                tftopo(ave_ersp,params.times,logfreqs,'limits', [params.times(1) params.times(end) logfreqs(1) logfreqs(end) -maxval maxval],...
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
                cbar;
                waitbar((k*n)/(len*Ncond),h_wait);
            end;
        end % Finish plotting all centroids for one condition
    end  % Finished all conditions
    if figureon
        delete(h_wait)
    end
end % Finished 'centroid' mode plot option

% cls_plotcompersp() - Commandline function, to visualizing cluster component ERSP images. 
%                    Displays the ERSP images of specified cluster components on separate figures,
%                    using one figure for all conditions. 
%                   The ERSPs can be visualized only if component ERSPs     
%                   were calculated and saved in the EEG datasets in the STUDY.
%                   These can be computed during pre-clustering using the GUI-based function
%                   pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                   and eeg_preclust(). A pop-function that calls this function is pop_clustedit().
% Usage:    
%                   >> [STUDY] = cls_plotcompersp(STUDY, ALLEEG, cluster, comps);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   cluster     - single cluster number.  
%
% Optional inputs:
%   comps      - [numeric vector]  -> indices of the cluster components to plot.
%                       'all'                       -> plot all the components in the cluster
%                                                      (as in cls_plotclustersp). {default: 'all'}.
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster ersp
%                     image mean, to allow quick replotting (unless cluster mean 
%                     already existed in the STUDY).  
%
%   Example:
%                         >> cluster = 4; comps= [1 7 10];  
%                         >> [STUDY] = cls_plotcompersp(STUDY,ALLEEG, cluster, comps);
%                    Plots components 1, 7 & 10  ersps of cluster 4 on separate figures. 
%
%  See also  pop_clustedit, pop_preclust, eeg_createdata, cls_plotclustersp         
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

% $ Log: cls_plotcompersp.m,v $

function STUDY = cls_plotcompersp(STUDY, ALLEEG, cls, varargin)
icadefs;
if ~exist('cls')
    error('cls_plotcompersp: you must provide a cluster number as an input.');
end
if isempty(cls)
   error('cls_plotcompersp: you must provide a cluster number as an input.');
end
if nargin == 3 % no components indices were given
    % Default plot all components of the cluster
    [STUDY] = cls_plotclustersp(STUDY, ALLEEG, 'clusters', cls, 'mode', 'comps');
    return
else
    comp_ind = varargin{1}; 
end
Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond = 1;
end
for ci = 1 : length(comp_ind) %for each comp
   rowcols(2) = ceil(sqrt(Ncond)); 
   rowcols(1) = ceil((Ncond)/rowcols(2));
   
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
        if ~isfield(ALLEEG(abset).etc,'icaerspparams')
            warndlg2([ 'Dataset ' num2str(abset) ' has no ERSP info, aborting'] , ['Abort - Plot ERSP']); 
            return;
        end
        params = ALLEEG(abset).etc.icaerspparams;
        sbplot(rowcols(1),rowcols(2),n), 
        if n == 1
            [ersp, logfreqs] = cls_readersp(ALLEEG, [STUDY.datasetinfo(STUDY.setind(:,STUDY.cluster(cls).sets(1,comp_ind(ci)))).index], comp);
            if isempty(ersp)
                warndlg2(['pop_clustedit: file '  ALLEEG(abset).etc.icalogersp ' was not found in path ' ALLEEG(abset).filepath], 'Abort - Plot ERSP' ); 
                return
            end
        end
        if Ncond >1
            a = [ 'ERSP, IC' num2str(comp) ' / ' subject ', ' STUDY.cluster(cls).name ', ' STUDY.condition{n} ];
        else
            a = ['ERSP, IC' num2str(comp) ' / ' subject  ', ' STUDY.cluster(cls(clus)).name];
        end
        tftopo(ersp(:,:,n),params.times,logfreqs,'limits', [params.times(1) params.times(end) logfreqs(1) logfreqs(end) -4 4],'title', a, 'verbose', 'off', 'axcopy', 'off');
        ft = str2num(get(gca,'yticklabel'));
        ft = exp(1).^ft;
        ft = unique(round(ft));
        ftick = get(gca,'ytick');
        ftick = exp(1).^ftick;
        ftick = unique(round(ftick));
        ftick = log(ftick);
        set(gca,'ytick',ftick);
        set(gca,'yticklabel', num2str(ft));
        cbar;
        axcopy(gcf, [' ft = str2num(get(gca,''''yticklabel'''')); ft = exp(1).^ft; ft = unique(round(ft)); fti = get(gca,''''ytick''''); fti = exp(1).^fti; fti = unique(round(fti));'...
            'fti = log(fti); set(gca, ''''ytick'''',fti); set(gca, ''''yticklabel'''',num2str(ft)); xlabel(''''Time [ms]''''); cbar; clear ft fti;' ]);
   end
end
