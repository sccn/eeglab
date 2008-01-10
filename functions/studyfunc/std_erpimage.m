% std_erpimage() - Plot an erpimage using multiple subject data
%
% Usage:
%   >> std_erpimage( STUDY, ALLEEG, typeplot, channel, projchan, title, ...
%                  smooth, decimate, sortingtype, sortingwin, ...
%                            sortingeventfield, renorm, options...);
% Inputs:
%   STUDY    - STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - global vector of EEG structures for the datasets included 
%              in the STUDY. ALLEEG for a STUDY set is typically created 
%              using load_ALLEEG().  
% Optional inputs:
%   'clusters' - [numeric vector|'all'] indices of clusters to plot.
%                If no component indices ('comps' below) are given, the average 
%                ERSPs of the requested clusters are plotted in the same figure, 
%                with ERSPs for different conditions (and groups if any) plotted 
%                in different colors. In 'comps' (below) mode, ERSP for each 
%                specified cluster are plotted in separate figures (one per 
%                condition), each overplotting cluster component ERSP plus the
%                average cluster ERSP in bold. Note this parameter has no effect 
%                if the 'comps' option (below) is used. {default: 'all'}
%   'comps'    - [numeric vector|'all'] indices of the cluster components to plot.
%                Note that 'comps', 'all' is equivalent to 'plotsubjects', 'on'.
%
% Optional inputs for channel plotting:
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel ERSP is plotted (using the 
%                same format as for the cluster component means described above)
%   'subject'  - [numeric vector]  In 'changrp' mode (above), index of 
%                the subject(s) to plot. Else by default, plot all components 
%                in the cluster.
%   'plotsubjects' - ['on'|'off'] When 'on', plot ERSP of all subjects.
%
% Commandline options:
%   'projchan' - Channel to back-project the selected component(s) to. 
%                If plotting channel activity, this argument is ignored. 
%                If [], the ICA component activation is plotted {default []}.
%   'title'      - ['string'] Plot title {default: []}
%   'smooth'     - Smoothing parameter (number of trials). {Default: 5} 
%                erpimage() equivalent: 'avewidth' 
%   'decimate'   - Decimate parameter (i.e. ratio of trials_in/trials_out).
%                erpaimge() equivalent: 'decimate' {Default: 0}
%   'sorttype'   - Sorting event type(s) ([int vector]; []=all). See Notes below.
%                Either a string or an integer.
%   'sortwin'    - Sorting event window [start, end] in seconds ([]=whole epoch)
%   'sortfield'  - Sorting field name. {default: none}. See Notes below.
%   options      - erpimage() options, separated by commas (Ex: 'erp', 'cbar'). 
%                {Default: none}. For further details see >> erpimage help   
%
% Outputs:
%  allphases   - all phase from erpimage
%  allerpimage - all sorting variables from erpimage
%  subjamp     - value of amplitude contribution [0 to 1] for each trial.
%                1 indicates that all subject contribute equally to the
%                amplitude of each trial. 0 indicates that only one
%                subject dominates all trials' amplitude.
%
% Author: Arnaud Delorme, SCCN & CERCO, CNRS, 2004-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Arnaud Delorme
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

%$Log: not supported by cvs2svn $
%Revision 1.4  2005/11/11 23:40:32  arno
%nothing
%
%Revision 1.3  2005/02/16 02:52:38  arno
%debug sumary array size
%
%Revision 1.2  2005/02/15 02:21:43  arno
%plot_indiv option
%
%Revision 1.1  2005/02/03 00:10:48  arno
%Initial revision
%

function [allphases, allsortvar, subjamptime, subjamptrial, globalent ] = std_erpimage( STUDY, ALLEEG, varargin);
    %subjind, data, sortvar, times, titleim, movewin, decim, varargin );
        
    if nargin < 3
        help std_erpimage;
        return;
    end;
    
    STUDY = pop_erpimparams(STUDY, 'default');

    [ opt moreparams ] = finputcheck( varargin, { ...
        'sorttype'    'string'  [] STUDY.etc.erpimparams.sortvar;
        'sortwin'     'string'  [] STUDY.etc.erpimparams.sortwin;
        'sortfield'   'string'  [] STUDY.etc.erpimparams.sortfield;
        'statistics'  'string'  [] STUDY.etc.erpimparams.statistics;
        'groupstats'  'string'  [] STUDY.etc.erpimparams.groupstats;
        'condstats'   'string'  [] STUDY.etc.erpimparams.condstats;
        'statmode'    'string'  [] STUDY.etc.erpimparams.statmode;
        'threshold'   'real'    [] STUDY.etc.erpimparams.threshold;
        'naccu'       'integer' [] STUDY.etc.erpimparams.naccu;
        'channels'    'cell'    []              {};
        'clusters'    'integer' []              [];
        'mode'        'string'  []              '';
        'comps'       {'integer','string'}  []              []; % for backward compatibility
        'plotsubjects' 'string' { 'on' 'off' }  'off';
        'subject'     'string'  []              '' }, ...
                                      'std_erpimage', 'ignore');
    if isstr(opt), error(opt); end;
    
    if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;

    if ~isempty(opt.subject), opt.groupstats = 'off'; disp('No group statistics for single subject'); end;
    if ~isempty(opt.subject), opt.condstats = 'off'; disp('No condition statistics for single subject'); end;
    
    if length(opt.comps) == 1
        opt.condstats = 'off'; opt.groupstats = 'off'; 
        disp('Statistics cannot be computed for single component');
    end;

    % read data from disk
    % -------------------
    if ~isempty(opt.channels)
         [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'data');
         STUDY = std_getepochevent(STUDY, ALLEEG, opt.sorttype, opt.sortwin, opt.softfield);
    else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', 'data');
    end;
    [STUDY.etc.datasort STUDY.etc.datind]  = std_getepochevent(STUDY, ALLEEG, opt.sorttype, opt.sortwin, opt.softfield);
    opt.legend = 'off';
    
    for index = 1:length(allinds)

        if length(allinds) > 1, try, subplot(nr,nc,index, 'align'); catch, subplot(nr,nc,index); end; end;
        if ~isempty(opt.channels)
            eval( [ 'alldata  = STUDY.changrp(allinds(index)).' opt.datatype 'data;' ]);
            eval( [ 'alltimes = STUDY.changrp(allinds(index)).' opt.datatype 'datatimes;' ]);
            setinds  = STUDY.changrp(allinds(index)).setinds;
        else
            eval( [ 'alldata  = STUDY.cluster(allinds(index)).' opt.datatype 'data;' ]);
            eval( [ 'alltimes = STUDY.cluster(allinds(index)).' opt.datatype 'datatimes;' ]);
            compinds = STUDY.cluster(allinds(index)).allinds;
            setinds  = STUDY.cluster(allinds(index)).setinds;
        end;
        
        % plot specific subject
        % ---------------------
        if ~isempty(opt.subject) & isempty(opt.comps)
            for c = 1:size(allersp,1)
                for g = 1:size(allersp,2)
                    for l=length(setinds{c,g}):-1:1
                        if ~strcmpi(opt.subject, STUDY.datasetinfo(setinds{c,g}(l)).subject)
                            allersp{c,g}(:,:,l) = [];
                        end;
                    end;
                end;
            end;
        end;
        
        % plot specific component
        % -----------------------
        if ~isempty(opt.comps)
            
            % find and select group
            % ---------------------
            sets   = STUDY.cluster(allinds(index)).sets(:,opt.comps);
            comps  = STUDY.cluster(allinds(index)).comps(opt.comps);
            grp    = STUDY.datasetinfo(sets(1)).group;
            grpind = strmatch( grp, STUDY.group );
            if isempty(grpind), grpind = 1; end;
            allersp = allersp(:,grpind);
            
            % find component
            % --------------
            for c = 1:size(allersp,1)
                for ind = length(compinds{1,grpind}):-1:1
                    if ~any(compinds{1,grpind}(ind) == comps) | ~any(setinds{1,grpind}(ind) == sets)
                        allersp{c}(:,:,ind) = [];
                    else
                        comp_names{c,1} = comps;
                    end;
                end;
            end;
            opt.subject = STUDY.datasetinfo(sets(1)).subject;
        end;
        
        % plot specific component
        % -----------------------
        if index == length(allinds), opt.legend = 'on'; end;
        erpimage(alldata, 
        
        
        [pgroup pcond pinter] = std_plottf(alltimes, allfreqs, allersp, 'condnames', STUDY.condition, 'subject', opt.subject, ...
                                           'legend', opt.legend, 'compinds', comp_names, 'datatype', opt.datatype,'plotmode', ...
                                           opt.plotmode, 'groupnames', STUDY.group, 'topovals', opt.plottf, 'unitx', 'Hz', ...
                                           'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
        if isempty(opt.channels), %title(sprintf('Cluster %d', allinds(index)));
            if length(allinds) > 1, 
                title([ STUDY.cluster(allinds(index)).name ' (' num2str(length(STUDY.cluster(allinds(index)).comps)),' ICs, ' ...
                        num2str(length(unique(STUDY.cluster(allinds(index)).sets(1,:)))) ' Ss)' ]);            
                set(gcf,'name','Cluster ERPIMAGE')
            elseif ~strcmp(opt.mode,'together') % if it is not the mean ERSP that is being shown (which is the case when 'cluster properties' is plotted then put cluster number on the corner of figure
                h = gca;
                axes('position',[0.04 0.96 0.1 0.06]); 
                text(0,0,[STUDY.cluster(allinds(index)).name],'fontsize',13 );
                axis off;
                if length(opt.comps) ~= 1
                    set(gcf,'name',['ERPIMAGE of ' STUDY.cluster(allinds(index)).name])
                else
                    set(gcf,'name',['ERPIMAGE of a Component from cluster ' STUDY.cluster(allinds(index)).name])
                end;
                axes(h);
            end;
        else                      
            title(sprintf('%s', opt.channels{index}));  
        end;
    end;
    
    
    
    plot_indiv = 1; % 0 or 1
    
    if nargin < 7
        help std_erpimage
        return;
    end;
    if iscell(movewin)
        erparg1 = movewin;
        erparg2 = decim;
        movewin = [ erparg1{1} erparg2{1} ];
        decim   = [ erparg1{2} erparg2{2} ];
    else
        erparg1 = [];
        erparg2 = [];
    end;
    
    % adapt parameters
    % ----------------    
    if size(data,1) == 1 & size(data,3) == 1
        data = reshape(data, size(data,2)/length(sortvar), length(sortvar));
    end;
    if length(movewin) == 1, movewin(2) = 1; end;
    if length(decim)   == 1, decim(2)   = 0;   end;
    
    % sorting
    % -------
    if nargin > 2 & ~isempty(sortvar)
        [tmpsort tmpind] = sort(sortvar);
        sortind          = subjind(tmpind);
    else
        sortind          = subjind;
    end;
    
    % plot all erpimage and retreive values
    % -------------------------------------
    values = unique(sortind);
    nsubj = length(values);
    realdecim = decim(1);
    if ~plot_indiv, extra_args = { 'noshow', 'yes' };
    else            extra_args = { };
    end;
    for ival = 1:nsubj
        idxsubj = find(subjind == values(ival));
        p(ival) = length(idxsubj);
        
        if plot_indiv, figure; end;
        if movewin(1) >= 0
            if isempty(erparg1)
                [ tmpdata tmpvar tmptrials tmplim axhndls,erp, ...
                  amps,cohers,cohsig,ampsig,outamps,phsangls,pshamp] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                              times, titleim, movewin(1), -decim(1), varargin{:}, extra_args{:});
                %{ times, titleim, movewin(1), -decim(1), varargin{:} }
            else
                [ tmpdata tmpvar tmptrials tmplim axhndls,erp, ...
                  amps,cohers,cohsig,ampsig,outamps,phsangls,pshamp] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                              times, titleim, movewin(1), -decim(1), erparg1{3:end}, extra_args{:});
                %{  times, titleim, movewin(1), -decim(1), erparg1{3:end} }
            end;
        else
            % subtract moving average from erpimage()
            % ---------------------------------------
            if isempty(erparg1)
                [ tmpdata tmpvar tmptrials tmplim axhndls,erp, ...
                  amps,cohers,cohsig,ampsig,outamps,phsangls,pshamp] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                                  times, titleim, -movewin(1), -decim(1), varargin{:}, extra_args{:});
                [ tmpdata2 tmpvar2 ] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                                  times, titleim, 1, 0, varargin{:}, 'noshow', 'yes');
            else
                [ tmpdata tmpvar tmptrials tmplim axhndls,erp, ...
                  amps,cohers,cohsig,ampsig,outamps,phsangls,pshamp] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                                  times, titleim, -movewin(1), -decim(1), erparg1{3:end}, extra_args{:});
                [ tmpdata2 tmpvar2 ] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                                  times, titleim, 1, 0, erparg1{3:end}, 'noshow', 'yes');
            end;
            if size(tmpdata2,2) > size(tmpdata,2)
                % randomly remove single trials
                % -----------------------------
                diflen   = size(tmpdata2,2) - size(tmpdata,2);
                indrm = find(shuffle([ ones(1, diflen) zeros(1, size(tmpdata,2)-diflen) ]));
                tmpdata2(:,indrm) = [];
                tmpvar2   (indrm) = [];
            end;
            maxlen = min(size(tmpdata,2), size(tmpdata2,2));
            tmpdata   = tmpdata2(:,1:maxlen)-tmpdata(:,1:maxlen);
            tmpvar    = tmpvar2(1:maxlen); % version non-smooth
            %tmpvar   = tmpvar(1:maxlen); % smooth version
            %phsangls  = phsangls(1:maxlen);
        end;

        % copy to variable regrouping all subject info
        % --------------------------------------------
        if ival == 1 
            realdecim  = size(tmpdata,2);
            outdata    = zeros(size(tmpdata,1), realdecim, nsubj);
            allsortvar = zeros(realdecim, nsubj);
            allphases  = zeros(realdecim, nsubj);
        end;
        
        if length(tmpdata) < size(outdata,2)
            outdata    = outdata   (:,1:length(tmpdata),:);
            allsortvar = allsortvar(1:length(tmpdata),:);
            if ~isempty(phsangls)
                allphases = allphases (1:length(tmpdata),:);
            end;
        end;
        minsize = min( size(outdata,2), size(tmpdata, 2) );
        outdata (:,1:minsize,ival) = tmpdata(:,1:minsize);     %/sum(sum(abs(tmpdata)));        
        allsortvar(1:minsize,ival) = tmpvar   (1:minsize)';
        
        if ~isempty(phsangls)
            tmp = movav(phsangls, [1:length(phsangls)], abs(movewin(1)), ...
                                    (length(phsangls)-abs(movewin(1))*2)/decim(1));
            allphases (:,ival) = tmp (1:size(allsortvar,1))';
        end;
    end;
    
    % uniform probability
    % -------------------
    uniform = p / sum(p);
    
    % resort trials and image them
    % ----------------------------
    newoutdata = mean(outdata,    3);
    newsortvar = mean(allsortvar, 2)';
    
    % plot global erpimage
    % --------------------
    timevect = linspace(tmplim(1), tmplim(2), size(outdata,1));
    if isempty(erparg2)
        % remove some arguments
        % ---------------------
        for index = length(varargin):-1:1
            if isstr(varargin{index}) & strcmpi(varargin{index}, 'align')
                varargin{index+1} = num2str(varargin{index+1});
            elseif isstr(varargin{index}) & strcmpi(varargin{index}, 'phasesort')
                varargin(index:index+1) = [];          
            elseif isstr(varargin{index}) & strcmpi(varargin{index}, 'plotamps')
                varargin(index) = [];          
            end;
        end;
        %{ times, titleim, movewin(2), -decim(2), varargin{:} }
        figure; erpimage( newoutdata, newsortvar, timevect, titleim, movewin(2), decim(2), 'nosort', varargin{:});
    else
        %{ times, titleim, movewin(2), -decim(2), erparg2{3:end} }
        figure; erpimage( newoutdata, newsortvar, timevect, titleim, movewin(2), decim(2), erparg2{3:end});
    end;
    return;
    
    % compute entropy
    % ---------------
    for sind = 1:nsubj % the loop below will select first
                       % the max trial among all subject
                       % second iteration will choose the second max...
                       % the loop has not been optimized
        
        % compute entropy array for this index
        % ------------------------------------
        subjcont = zeros(1, length(subjind));
        arrayp   = zeros(size(outdata,1), realdecim, nsubj);
        for index = 1:realdecim
            
            trials = outdata(:, index:realdecim:end); % nsubj trials (1 per subject)
            
            for indtime = 1:size(trials,1)
                [tmp indorder] = sort(trials(indtime,:));
                arrayp(indtime, index, indorder(sind)) = arrayp(indtime, index, indorder(1))+1;
            end;
            
        end;
        
        % compute entropy across time
        % ---------------------------
        arrayptime = squeeze(sum(arrayp,2));
        for indtime = 1:size(trials,1)
            p = arrayptime(indtime,:);
            p = p / sum(p);
            p(find(p == 0)) = [];
            subjamptime(indtime,sind) = -sum(p.*log(p)) / -sum(uniform.*log(uniform)) ;
        end;
        
        % compute entropy across trials
        % -----------------------------
        arrayptrial = squeeze(sum(arrayp,1));
        for indtrial = 1:size(arrayptrial,1)
            p = arrayptime(indtrial,:);
            p = p / sum(p);
            p(find(p == 0)) = [];
            subjamptrial(indtrial,sind) = -sum(p.*log(p)) / -sum(uniform.*log(uniform)) ;
        end;
    
        % global entropy
        % --------------
        p = squeeze(sum(sum(arrayp,1),2));
        p = p / sum(p);
        p(find(p == 0)) = [];
        globalent(sind) =  -sum(p.*log(p)) / -sum(uniform.*log(uniform));
    end;
    
    figure; plot(timevect, mean(subjamptime,2));
    figure; plot(1:realdecim, mean(subjamptrial,2));
    
    return
    
