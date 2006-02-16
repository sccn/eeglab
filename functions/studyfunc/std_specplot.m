% cls_plotclustspec() - Commandline function, to visualizing cluster/s components spectra. 
%                   Either displays mean spectra of all requested clusters in the same figure, 
%                   with spectra for different conditions (if any) plotted in different colors. 
%                   Or displays spectra for each specified cluster in separate figures (per condition),  
%                   each containing the cluster components plus the average cluster spectrum in bold.
%                   The spectra can be visualized only if component spectra     
%                   were calculated and saved in the EEG datasets in the STUDY.
%                   These can be computed during pre-clustering using the GUI-based function
%                   pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                   and eeg_preclust(). A pop-function that calls this function is pop_clustedit().
% Usage:    
%                   >> [STUDY] = cls_plotclustspec(STUDY, ALLEEG, key1, val1, key2, val2);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Optional inputs:
%   'clusters'   - [numeric vector]  -> specific cluster numbers to plot.
%                     'all'                         -> plot all clusters in STUDY.
%                     {default: 'all'}.
%   'mode'       - ['centroid'|'comps'] a plotting mode. In 'centroid' mode, the average spectra 
%                     of the requested clusters are plotted in the same figure, with spectra for  
%                     different conditions (if any) plotted in different colors. In 'comps' mode, spectra
%                     for each specified cluster are plotted in separate figures (per condition), each 
%                     containing  cluster components spectra plus the average cluster spectrum in bold.
%                     {default: 'centroid'}.
%   'figure'       - ['on'|'off'] for the 'centroid' mode option, plots on
%                     a new figure ('on')  or plots on current figure ('off').
%                     {default: 'on'}.
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster 
%                     mean spectra, to allow quick replotting (unless cluster means 
%                     already exists in the STUDY).  
%
%   Example:
%                         >> [STUDY] = cls_plotclustspec(STUDY,ALLEEG, 'clusters', 2, 'mode', 'comps');
%                    Plots cluster 2 components spectra along with the mean spectra in bold. 
%
%  See also  pop_clustedit, pop_preclust, eeg_createdata, cls_plotcompspec         
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
% Revision 1.3  2006/02/16 00:01:21  arno
% fixing axis limits
%
% Revision 1.2  2006/02/15 22:59:20  arno
% adding scaling etc...
%

function STUDY = cls_plotclustspec(STUDY, ALLEEG,  varargin)
icadefs;
% Set default values
cls = 2:length(STUDY.cluster); % plot all clusters in STUDY
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
                    error('cls_plotclustersp: ''clusters'' input takes either specific clusters (numeric vector) or keyword ''all''.');
                end
            end
        case 'mode' % Plotting mode 'centroid' / 'comps'
            mode = varargin{k-1};
        case 'figure'
            if strcmpi(varargin{k-1},'off')
                figureon = 0;
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
    Ncond =1;
end
% Plot all the components in the cluster ('comps' mode)
if strcmpi(mode, 'comps')         
    for clus = 1: length(cls) % For each cluster requested
        len = length(STUDY.cluster(cls(clus)).comps);
        if ~isfield(STUDY.cluster(cls(clus)).centroid,'spec')
            STUDY = cls_centroid(STUDY,ALLEEG, cls(clus) , 'spec');
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
                clusnval = cls_clusread(STUDY, ALLEEG, cls(clus),'spec',n);
            catch,
                warndlg2([ 'Some spectra information is missing, aborting'] , ['Abort - Plot spectra'] );   
                return;
            end
            handl(n) = sbplot(rowcols(1),rowcols(2),n);
            ave_spec = STUDY.cluster(cls(clus)).centroid.spec{n};
            f = STUDY.cluster(cls(clus)).centroid.spec_f;
            plot(f,clusnval.spec,'color', [0.5 0.5 0.5]);
            hold on
            plot(f,ave_spec,'k','linewidth',2);
            axis([f(1) f(end) -10 20]);
            xlabel('Frequency [Hz]');
            ylabel('Power [dB]');
            title(['Spectra, '  STUDY.cluster(cls(clus)).name ', ' STUDY.condition{n} ', ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) 'Ss']);
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
                    axis([f(1) f(end)  ylimits(1)  ylimits(2) ]);
                    axcopy;
                end
            end
        end % finished one condition
    end % finished all requested clusters 
end % Finished 'comps' mode plot option
       
% Plot clusters mean spec
if strcmpi(mode, 'centroid') 
    len = length(cls);
    rowcols(2) = ceil(sqrt(len)); rowcols(1) = ceil((len)/rowcols(2));
    if figureon
        try 
            % optional 'CreateCancelBtn', 'delete(gcbf); error(''USER ABORT'');', 
            h_wait = waitbar(0,['Computing Spectra ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,['Computing Spectra ...'],'position', [300, 200, 300, 48]);
        end
        figure
	end
    color_codes = {'b', 'r', 'g', 'c', 'm', 'y', 'k','b--', 'r--', 'g--', 'c--', 'm--', 'y--', 'k--','b-.', 'r-.', 'g-.', 'c-.', 'm-.', 'y-.', 'k-.'};
    orient tall
    
    min_spec = Inf;
    max_spec = -Inf;
    for k = 1:len % Go through the clusters
        if ~isfield(STUDY.cluster(cls(k)).centroid,'spec')
            STUDY = cls_centroid(STUDY,ALLEEG, cls(k) , 'spec');
        end
        if len ~= 1
            sbplot(rowcols(1),rowcols(2),k) ; 
        end
        hold on;
        for n = 1:Ncond
            if k == 1
                leg_color{n} = [STUDY.condition{n}];
            end
            ave_spec = STUDY.cluster(cls(k)).centroid.spec{n};
            f = STUDY.cluster(cls(k)).centroid.spec_f;
            plot(f,ave_spec,color_codes{n},'linewidth',2);
            
            min_spec = min(min(ave_spec), min_spec);
            max_spec = max(max(ave_spec), max_spec);
            
            if n == Ncond
                a = [ STUDY.cluster(cls(k)).name ', '  num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ];
                title(a);
                set(gcf,'Color', BACKCOLOR);
                set(gca,'UserData', leg_color);
                set(gcf,'UserData', leg_color);
                axis([f(1) f(end) -10 20]);
                axcopy(gcf, 'leg_color = get(gca,''''UserData'''') ; legend(leg_color); xlabel(''''Frequency [Hz]'''');ylabel(''''Power [dB]'''') ;');
                if figureon
                    waitbar(k/len,h_wait);
                end
            end
            if (k == len) & (n == Ncond)
                xlabel('Frequency [Hz]');
                ylabel('Power [dB]');
                if len ~= 1
                    maintitle = ['Average spectra for several clusters across all conditions'];
                    a = textsc(maintitle, 'title'); 
                    set(a, 'fontweight', 'bold'); 
                else
                    a = [ STUDY.cluster(cls(k)).name ' spectra, '  num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ];
                    title(a);
                end
                diff_spec = max_spec-min_spec;
                ylim( [ min_spec - 0.1*diff_spec max_spec + 0.1*diff_spec]);
                set(gcf,'Color', BACKCOLOR);
                legend(leg_color);
                if figureon
                    delete(h_wait);
                end
            end
        end % finished the different conditions
    end % finished all clusters 
end % finished 'centroid' plot mode

