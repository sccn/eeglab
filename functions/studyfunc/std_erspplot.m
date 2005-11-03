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
%                     'all'                         -> plot all clusters in STUDY.
%                     {default: 'all'}.
%   'mode'       - ['centroid'|'comps'] a plotting mode. In 'centroid' mode, the average ERSPs 
%                     of the requested clusters are plotted in the same figure - one per condition. 
%                     In 'comps' mode, component ERSPs for each cluster are plotted in a
%                     separate figure (per condition) with the cluster mean ERSP.
%                     {default: 'centroid'}.
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
            tftopo(ave_ersp,ersp_times,logfreqs,'limits', [ersp_times(1) ersp_times(end) logfreqs(1) logfreqs(end) -lim lim],...
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
        rowcols(2) = ceil(sqrt(len)); rowcols(1) = ceil((len)/rowcols(2));
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
    lim(1:Ncond) = 0;   
    % Compute cluster centroid
    for k = 1:len 
        if ~isfield(STUDY.cluster(cls(k)).centroid,'ersp')
            STUDY = cls_centroid(STUDY,ALLEEG, cls(k) , 'ersp');
        end
        if isempty(STUDY.cluster(cls(k)).centroid.ersp)
            warndlg2(['eeg_clustedit: some .icaersp files could not be found for cluster ' STUDY.cluster(cls(k)).name ], 'Abort - Plot ERSP' ); 
            return
        end
        for n = 1:Ncond
            lim(n) = lim(n)+STUDY.cluster(cls(k)).centroid.ersp_limits{n};%plotting limits
        end
    end
    lim = lim./len;
    params = ALLEEG(STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(cls(1)).sets(1,1))).index).etc.icaerspparams;            
    for n = 1:Ncond
        if figureon 
            figure
            maintitle = ['Average ERSP for all clusters, ' STUDY.condition{n}];
            a = textsc(maintitle, 'title'); 
            set(a, 'fontweight', 'bold'); 
        end
        orient tall
        set(gcf,'Color', BACKCOLOR);
        for k = 1:len 
            if len ~= 1
                sbplot(rowcols(1),rowcols(2),k) , 
            end
            a = [ STUDY.cluster(cls(k)).name ' ERSP, ' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss, ' STUDY.condition{n}];
            ave_ersp = STUDY.cluster(cls(k)).centroid.ersp{n};
            logfreqs = STUDY.cluster(cls(k)).centroid.ersp_logf;
            if figureon % plot on a new figure                
                tftopo(ave_ersp,params.times,logfreqs,'limits', [params.times(1) params.times(end) logfreqs(1) logfreqs(end) -lim(n) lim(n)],...
                    'title', a, 'verbose', 'off', 'axcopy', 'off');
                axcopy(gcf, [' ft = str2num(get(gca,''''yticklabel'''')); ft = exp(1).^ft; ft = unique(round(ft)); fti = get(gca,''''ytick''''); fti = exp(1).^fti; fti = unique(round(fti));'...
                     'fti = log(fti); set(gca, ''''ytick'''',fti); set(gca, ''''yticklabel'''',num2str(ft)); xlabel(''''Time [ms]''''); ylabel(''''Frequency [Hz]''''); cbar; clear ft fti;' ]);
                if k ~= len
                    set(gca, 'xtick', [], 'ytick', []);
                    xlabel('');
                    ylabel('');
                else
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
                end
                waitbar((k*n)/(len*Ncond),h_wait);
            else % plot on existing figure (plot only one condition)
               ersp{n}.ersp = ave_ersp;
               ersp{n}.times = params.times;
               ersp{n}.logf = logfreqs;
               ersp{n}.limits = [params.times(1) params.times(end) logfreqs(1) logfreqs(end) -mean(lim) mean(lim)];
               ersp{n}.title = a;
               if n == Ncond % finished computing all conditions
                   tftopo(ersp{1}.ersp,ersp{1}.times,ersp{1}.logf,'limits', ersp{1}.limits, 'title', ersp{1}.title, 'verbose', 'off', 'axcopy', 'off');
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
                   set(get(gca,'Title'),'FontSize',10);
                   set(get(gca,'Ylabel'),'FontSize',10);
                   set(get(gca,'Xlabel'),'FontSize',10);   
                   set(get(gcf,'Children'),'FontSize',10);
                   set(gcf, 'UserData', ersp);
                   set(gca,'UserData', ersp);
                   axcopy(gcf, ['ersp = get(gca, ''''UserData''''); len = length(ersp); rowcols(2) = ceil(sqrt(len)); rowcols(1) = ceil((len)/rowcols(2)); '...
                   ' for k = 1:len, subplot(rowcols(1),rowcols(2),k), '...
                   ' tftopo(ersp{k}.ersp,ersp{k}.times,ersp{k}.logf,''''limits'''', ersp{k}.limits, ''''title'''', ersp{k}.title, ''''verbose'''', ''''off'''', ''''axcopy'''', ''''off'''');' ...
                   ' ft = str2num(get(gca,''''yticklabel'''')); ft = exp(1).^ft; ft = unique(round(ft)); fti = get(gca,''''ytick''''); fti = exp(1).^fti; fti = unique(round(fti));'...
                   'fti = log(fti); set(gca, ''''ytick'''',fti); set(gca, ''''yticklabel'''',num2str(ft)); xlabel(''''Time [ms]'''');' ...
                   'ylabel(''''Frequency [Hz]''''); end; cbar; clear ft fti ersp;' ]);
               end
            end                
        end % Finish plotting all centroids for one condition
    end  % Finished all conditions
    if figureon
        delete(h_wait)
    end
end % Finished 'centroid' mode plot option
