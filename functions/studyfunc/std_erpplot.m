% std_erpplot() - Command line function to plot STUDY cluster component ERPs. Either 
%                 displays grand mean ERPs for all requested clusters in the same figure, 
%                 with ERPs for different conditions (if any) plotted in different colors. 
%                 Else, displays ERP for each specified cluster in separate figures 
%                 (per condition), each containing the cluster component ERPs plus 
%                 the grand mean cluster ERP (in bold). ERPs can be plotted only if 
%                 component ERPs were computed and saved in the STUDY EEG
%                 datasets. 
%                 These can be computed during pre-clustering using the gui-based 
%                 function pop_preclust() or the equivalent command line functions 
%                 eeg_createdata() and eeg_preclust(). Called by pop_clustedit().
%                 and std_propplot().
% Usage:    
%   >> [STUDY] = std_erpplot(STUDY, ALLEEG, key1, val1, key2, val2);  
%   >> [STUDY erpdata erptimes pgroup pcond pinter] = std_erpplot(STUDY, ALLEEG, ...);  
%
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the datasets included 
%                in the STUDY. A STUDY set ALLEEG is typically created by load_ALLEEG().  
% Optional inputs for channel plotting:
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel ERP is plotted (using the 
%                same format as for the cluster component means described
%                above). Default is to plot all channels.
%   'subject'  - [numeric vector]  In 'changrp' mode (above), index of 
%                the subject(s) to plot. Else by default, plot all components 
%                in the cluster.
%   'plotsubjects' - ['on'|'off'] When 'on', plot ERP of all subjects.
%   'noplot'   - ['on'|'off'] When 'on', only return output values. Default
%                is 'off'.
%   'topoplotopt' - [cell array] options for topoplot plotting.
%
% Optional inputs for component plotting:
%   'clusters' - [numeric vector|'all'] indices of clusters to plot.
%                If no component indices ('comps' below) are given, the average 
%                ERPs of the requested clusters are plotted in the same figure, 
%                with ERPs for different conditions (and groups if any) plotted 
%                in different colors. In 'comps' (below) mode, ERPS for each 
%                specified cluster are plotted in separate figures (one per 
%                condition), each overplotting cluster component ERPs plus the
%                average cluster ERP in bold. Note this parameter has no effect 
%                if the 'comps' option (below) is used. {default: 'all'}
%   'comps'    - [numeric vector|'all'] indices of the cluster components to plot.
%                Note that 'comps', 'all' is equivalent to 'plotsubjects', 'on'.
%
% Other optional inputs:
%   'key','val' - All optional inputs to pop_erpparams() are also accepted here
%                 to plot subset of time, statistics etc. The values used by default
%                 are the ones set using pop_erpparams() and stored in the
%                 STUDY structure.
%
% Outputs:
%   STUDY      - the input STUDY set structure with plotted cluster mean
%                ERPs data to allow quick replotting 
%   erpdata    - [cell] ERP data for each condition, group and subjects.
%                size of cell array is [nconds x ngroups]. Size of each element
%                is [times x subjects] for data channels or [times x components]
%                for component clusters. This array may be gicen as input 
%                directly to the statcond() function or std_stats function
%                to compute statistics.
%   erptimes   - [array] ERP time point latencies.
%   pgroup     - [array or cell] p-values group statistics. Output of the 
%                statcond() function.
%   pcond      - [array or cell] condition statistics. Output of the statcond() 
%                function.
%   pinter     - [array or cell] groups x conditions statistics. Output of
%                statcond() function.
%
%   Example:
%            >> [STUDY] = std_erpplot(STUDY,ALLEEG, 'clusters', 2, 'comps', 'all');
%               % Plot cluster-2 component ERPs plus the mean ERP in bold. 
%
%  See also  pop_clustedit(), pop_preclust(), eeg_createdata(), eeg_preclust(). std_propplot()
%
% Authors: Arnaud Delorme, CERCO, August, 2006-

% Copyright (C) Arnaud Delorme, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [STUDY, erpdata, alltimes, pgroup, pcond, pinter] = std_erpplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_erpplot;
    return;
end
erpdata = []; alltimes = []; 
pgroup = []; pcond = []; pinter = [];

% find datatype and default options
% ---------------------------------
dtype = 'erp';
dsubtype = '';
for ind = 1:2:length(varargin)
    if strcmpi(varargin{ind}, 'datatype')
        dtype = varargin{ind+1}; 
    end
end
if strcmpi(dtype(1:3), 'erp' )
    if length(dtype) > 3, dsubtype = dtype(4:end); dtype = 'erp'; end
elseif strcmpi(dtype(1:4), 'spec')
    if length(dtype) > 4, dsubtype = dtype(5:end); dtype = 'spec'; end
end

% get parameters
% --------------
eval( [ 'tmp = pop_' dtype 'params(STUDY, varargin{:});' ...
        'params = tmp.etc.' dtype 'params; clear tmp;' ] );
statstruct.etc = STUDY.etc; 
statstruct = pop_statparams(statstruct, varargin{:});

% potentially missing fields
% --------------------------
fields     = { 'filter' 'subtractsubjectmean' 'timerange' 'freqrange' 'topotime' 'topofreq' 'averagechan'};
defaultval = { [] 'off' [] [] [] [] };
for ind=1:length(fields)
    if ~isfield(params, fields{ind})
        params = setfield(params, fields{ind}, defaultval{ind}); 
    end
end

% decode parameters
% -----------------
if isempty(varargin)
    tmplocs = eeg_mergelocs(ALLEEG.chanlocs);
    options.channels = { tmplocs.labels };
else
    options = mystruct(varargin);
end
options = myrmfield( options, myfieldnames(params));
options = myrmfield( options, myfieldnames(statstruct.etc.statistics));
options = myrmfield( options, { 'threshold' 'statistics' } ); % for backward compatibility
opt = finputcheck( options, ...
                             { 'design'      'integer' []              STUDY.currentdesign;
                               'plotstderr'  'string'  []              'off';
                               'channels'    'cell'    []              {};
                               'clusters'    'integer' []              [];
                               'datatype'    'string'  {}              'erp';
                               'mode'        'string'  []              ''; % for backward compatibility (now used for statistics)
                               'comps'       { 'string','integer' } [] []; % for backward compatibility
                               'statmode'    'string'  { 'subjects','common','trials' } 'subjects'; % ignored
                               'plotmode'    'string'  { 'normal','condensed' }  'normal';
                               'unitx'       'string'  { 'ms','Hz' }    'ms';
                               'plotsubjects' 'string' { 'on','off' }  'off';
                               'detachplots' 'string'  { 'on','off' }  params.detachplots;
                               'noplot'      'string'  { 'on','off' }  'off';
                               'topoplotopt' 'cell'    {}              { 'style' 'both' };
                               'subject'     'string'  []              '' }, 'std_erpplot');
if ischar(opt), error(opt); end
if ischar(opt.comps), opt.comps = []; opt.plotsubjects = 'on'; end
if ~isempty(params.topofreq) && strcmpi(opt.datatype, 'spec'),  params.topotime  = params.topofreq; end
if ~isempty(params.freqrange), params.timerange = params.freqrange; end
datatypestr = upper(opt.datatype);
if strcmpi(datatypestr, 'spec'), datatypestr = 'Spectrum'; end
    
% =======================================================================
% below this line, all the code should be non-specific to ERP or spectrum
% =======================================================================
allconditions  = {};
allgroups      = {};
condname       = '';
groupname      = '';
if length(STUDY.design(opt.design).variable) > 0, allconditions = STUDY.design(opt.design).variable(1).value; condname  = STUDY.design(opt.design).variable(1).label; end
if length(STUDY.design(opt.design).variable) > 1, allgroups     = STUDY.design(opt.design).variable(2).value; groupname = STUDY.design(opt.design).variable(2).label;  end

% for backward compatibility
% --------------------------
stats = statstruct.etc.statistics;
stats.fieldtrip.channelneighbor = struct([]); % assumes one channel or 1 component
if isempty(STUDY.design(opt.design).variable)
    stats.paired = { };
else
    stats.paired = { STUDY.design(opt.design).variable(:).pairing };
end
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end
if strcmpi(stats.singletrials, 'off') && ((~isempty(opt.subject) || ~isempty(opt.comps)))
    if strcmpi(stats.condstats, 'on') || strcmpi(stats.groupstats, 'on')
        stats.groupstats = 'off';
        stats.condstats   = 'off'; 
        disp('No statistics for single subject/component, to get statistics compute single-trial measures'); 
    end
end

if ~isempty(params.topotime) && ~isnan(params.topotime(1)) && length(opt.channels) < 5 && isempty(opt.clusters)
    warndlg2(strvcat('ERP parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
    return;
end

plotcurveopt = {};
if length(opt.clusters) > 1
    plotcurveopt = { 'figure' 'off' }; 
    params.plotconditions = 'together';
    params.plotgroups     = 'together';
    stats.condstats  = 'off'; 
    stats.groupstats = 'off';
end
% if length(opt.channels) > 1 && strcmpi(opt.plotconditions, 'together') && strcmpi(opt.plotgroups, 'together')
%     plotcurveopt = { 'figure' 'off' }; 
%     opt.plotconditions = 'together';
%     opt.plotgroups     = 'together';
%     opt.condstats  = 'off'; 
%     opt.groupstats = 'off';
% end
alpha    = fastif(strcmpi(stats.mode, 'eeglab'), stats.eeglab.alpha, stats.fieldtrip.alpha);
mcorrect = fastif(strcmpi(stats.mode, 'eeglab'), stats.eeglab.mcorrect, stats.fieldtrip.mcorrect);
method   = fastif(strcmpi(stats.mode, 'eeglab'), stats.eeglab.method, ['Fieldtrip ' stats.fieldtrip.method ]);
plotcurveopt = { plotcurveopt{:} ...
   'ylim',           params.ylim, ...
   'threshold',      alpha ...
   'unitx'           opt.unitx, ...
   'filter',         params.filter, ...
   'plotgroups',     params.plotgroups, ...
   'effect',         stats.effect, ...
   'plotconditions', params.plotconditions };

% channel plotting
% ----------------
axcopyflag = 1;
if ~isempty(opt.channels)
    if (isempty(params.topotime) || isnan(params.topotime(1))) && length(opt.channels) > 1 && strcmpi(stats.singletrials, 'on')
        error('Cannot plot several channels on the same figure when using single trial statistics');
    end

    chaninds = 1:length(opt.channels);

    if strcmpi(opt.datatype, 'erp')
        [STUDY, erpdata, alltimes, ~, ~, fileparams] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels(chaninds), 'timerange', params.timerange, ...
                'subject', opt.subject, 'singletrials', stats.singletrials, 'design', opt.design, 'datatype', [dtype dsubtype]);
    else
        [STUDY, erpdata, alltimes, ~, ~, fileparams] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels(chaninds), 'freqrange', params.freqrange, ...
                'subject', opt.subject, 'singletrials', stats.singletrials, 'design', opt.design, 'datatype', [dtype dsubtype], 'rmsubjmean', params.subtractsubjectmean);
    end
    if isfield(fileparams, 'specmode') && ~strcmpi(fileparams.specmode, 'fft'), opt.unitx = [ opt.unitx 'psd' ]; end
    if ~strcmpi(params.averagechan, 'off') && length(chaninds) > 1
        for index = 1:length(erpdata(:))
            if strcmpi(params.averagechan, 'on')
                erpdata{index} = squeeze(mean(erpdata{index},2));
            else
                erpdata{index} = squeeze(sqrt(mean(erpdata{index}.^2,2)));
                opt.unitx = 'rmsms';
            end
        end
        if strcmpi(params.averagechan, 'rms')
            opt.unitx = [ 'rms' opt.unitx ];
        end
    end
    if isempty(erpdata), return; end

    % select specific time
    % --------------------
    if ~isempty(params.topotime) && ~isnan(params.topotime(1))
        [tmp, ti1] = min(abs(alltimes-params.topotime(1)));
        [tmp, ti2] = min(abs(alltimes-params.topotime(end)));
        for condind = 1:length(erpdata(:))
            if ~isempty(erpdata{condind})
                erpdata{condind} = mean(erpdata{condind}(ti1:ti2,:,:),1);
            end
        end
    end

    % compute baseline for spectrum
    % -----------------------------
    if strcmpi(params.subtractsubjectmean, 'on') && strcmpi(opt.datatype, 'spec')
        count    = 0;
        for iSpec = 1:length(erpdata(:))
            if ~isempty(erpdata{iSpec})
                if count == 0, meanspec = zeros(size(erpdata{iSpec},1),1); end
                if strcmpi(stats.singletrials, 'on')
                    meanspec = meanspec + mean(erpdata{iSpec},2);
                else
                    meanspec = meanspec + erpdata{iSpec};
                end
                count = count+1;
            end
        end
        meanspec = meanspec/count;
        erpdata  = cellfun(@(x)bsxfun(@minus, x, meanspec), erpdata, 'uniformoutput', false);
    end

    % compute statistics
    % ------------------
    if (isempty(params.topotime) || any(isnan(params.topotime))) && length(alpha) > 1
        alpha = alpha(1);
    end
    if ~isempty(params.topotime) && all(~isnan(params.topotime))
         statstruct = std_prepare_neighbors(statstruct, ALLEEG, 'channels', opt.channels);
         stats.fieldtrip.channelneighbor = statstruct.etc.statistics.fieldtrip.channelneighbor;
    end
    [pcond, pgroup, pinter] = std_stat(erpdata, stats);
    if (~isempty(pcond) && length(pcond{1}) == 1) || (~isempty(pgroup) && length(pgroup{1}) == 1), pcond = {}; pgroup = {}; pinter = {}; end % single subject STUDY                                
    if ~isempty(params.topotime) && length(opt.channels) > 5 && ndims(erpdata{1}) < 3, pcond = {}; pgroup = {}; pinter = {}; end % topo plotting for single subject
    if strcmpi(opt.noplot, 'on') return; end
    
    % get titles (not included in std_erspplot because it is not possible
    % to merge channels for that function
    % -----------------------------------
    locsOri = eeg_mergelocs(ALLEEG.chanlocs);
    locs = locsOri(std_chaninds(STUDY, opt.channels(chaninds)));
    if ~strcmpi(params.averagechan, 'off') && length(chaninds) > 1
        if length(chaninds) ~= length(locsOri)
            chanlabels = { locs.labels };
            chanlabels(2,:) = {','};
            chanlabels(2,end) = {''};
            locs(1).labels = [ chanlabels{:} ];
        else
            locs(1).labels = 'All channels';
        end
        locs(2:end) = [];    
    end
    [alltitles, alllegends ] = std_figtitle('threshold', alpha, 'mcorrect', mcorrect, 'condstat', stats.condstats, 'cond2stat', stats.groupstats, ...
                             'statistics', method, 'condnames', allconditions, 'plotsubjects', opt.plotsubjects, 'cond2names', allgroups, 'chanlabels', { locs.labels }, ...
                             'subject', opt.subject, 'valsunit', opt.unitx, 'vals', params.topotime, 'datatype', datatypestr, 'cond2group', params.plotgroups, ...
                             'condgroup', params.plotconditions, 'effect', stats.effect, 'factor1', condname, 'factor2', groupname);

    % plot
    % ----
    indNonEmpty = find(~cellfun(@isempty, erpdata(:)));
    if ~isreal(erpdata{indNonEmpty(1)}(1)) % for spectrum FFT data
        tmperpdata = cellfun(@(x)x.*conj(x), erpdata, 'uniformoutput', false);
    else
        tmperpdata = erpdata;
    end
    if ~isempty(params.topotime) && all(~isnan(params.topotime))
        std_chantopo(tmperpdata, 'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'caxis', params.ylim, ...
                                      'chanlocs', locs, 'threshold', alpha, 'titles', alltitles, 'topoplotopt', opt.topoplotopt, 'effect', stats.effect);
    else
        std_plotcurve(alltimes, tmperpdata, 'groupstats', pgroup, 'legend', alllegends, 'condstats', pcond, 'interstats', pinter, ...
            'chanlocs', locs, 'titles', alltitles, 'plotsubjects', opt.plotsubjects, 'plotstderr', opt.plotstderr, ...
            'condnames', allconditions, 'groupnames', allgroups, plotcurveopt{:});
    end

    set(gcf,'name',['Channel ' datatypestr ]);
    axcopy(gca);
else 
    if length(opt.clusters) > 1 && strcmpi(stats.singletrials, 'on')
        error('Cannot plot several components on the same figure when using single trial statistics');
    end
    
   % plot component
    % --------------
    if length(opt.clusters) > 1, figure('color', 'w'); end
    nc = ceil(sqrt(length(opt.clusters)));
    nr = ceil(length(opt.clusters)/nc);
    comp_names = {};
    
    for index = 1:length(opt.clusters)

        if length(opt.clusters) > 1, subplot(nr,nc,index); end
        if strcmpi(opt.datatype, 'erp')
            [STUDY, erpdata, alltimes, ~, ~, fileparams] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters(index), 'timerange', params.timerange, ...
                    'component', opt.comps, 'singletrials', stats.singletrials, 'design', opt.design, 'datatype', [dtype dsubtype]);
        else
            [STUDY, erpdata, alltimes, ~, ~, fileparams] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters(index), 'freqrange', params.freqrange, ...
                    'component', opt.comps, 'singletrials', stats.singletrials, 'design', opt.design, 'datatype', [dtype dsubtype], 'rmsubjmean', params.subtractsubjectmean);
        end
        if isfield(fileparams, 'specmode') && ~strcmpi(fileparams.specmode, 'fft'), opt.unitx = [ opt.unitx 'psd' ]; end
        if isempty(erpdata), return; end

        % plot specific component
        % -----------------------
        if ~isempty(opt.comps)
            comp_names = { STUDY.cluster(opt.clusters(index)).comps(opt.comps) };
            opt.subject = STUDY.datasetinfo(STUDY.cluster(opt.clusters(index)).sets(1,opt.comps)).subject;
        end
        
        % remove NaNs and generate labels
        % -------------------------------
        erpdata2 = erpdata;
        subjects       = { STUDY.datasetinfo(STUDY.cluster(opt.clusters(index)).sets(1,:)).subject };
        for iDat = 1:length(erpdata2(:))
            if all(cellfun(@(x)size(x,2), erpdata2(:)) == length(subjects))  % NOT single trial data
                keepInd = find(~isnan(erpdata2{iDat}(1,:)));
                erpdata2{iDat} = erpdata2{iDat}(:,keepInd);
                tmpSubjects = subjects(keepInd);
                comps       = STUDY.cluster(opt.clusters(index)).comps(keepInd);
            else
                tmpSubjects = subjects;
                comps       = STUDY.cluster(opt.clusters(index)).comps;
            end
            for iKeep = 1:length(tmpSubjects)
                sbtitles{iDat}{iKeep} = [ tmpSubjects{iKeep}  '/IC' num2str(comps(iKeep)) ];
            end
        end
        sbtitles = reshape(sbtitles, size(erpdata2));
        
        [pcond, pgroup, pinter] = std_stat(erpdata2, stats);
        if strcmpi(opt.noplot, 'on'), return; end
            
        [alltitles, alllegends ] = std_figtitle('threshold', alpha, 'plotsubjects', opt.plotsubjects, 'mcorrect', mcorrect, 'condstat', stats.condstats, 'cond2stat', stats.groupstats, ...
                                 'statistics', method, 'condnames', allconditions, 'cond2names', allgroups, 'clustname', STUDY.cluster(opt.clusters(index)).name, 'compnames', comp_names, ...
                                 'subject', opt.subject, 'valsunit', opt.unitx, 'vals', params.topotime, 'datatype', datatypestr, 'cond2group', params.plotgroups, 'condgroup', params.plotconditions, ...
                                 'effect', stats.effect, 'factor1', condname, 'factor2', groupname);
        
        if length(opt.clusters) > 1 && index < length(opt.clusters), alllegends = {}; end
        std_plotcurve(alltimes, erpdata2, 'condnames', allconditions, 'legend', alllegends, 'groupnames', allgroups, 'plotstderr', opt.plotstderr, ...
                                          'titles', alltitles, 'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                          'plotsubjects', opt.plotsubjects, plotcurveopt{:});
%--------------------------------------------------------------------------                                      
        if all([strcmp(opt.plotsubjects,'on') strcmp(opt.detachplots,'on')])
            
             [alltitlestmp] = std_figtitle('threshold', alpha, 'plotsubjects', opt.plotsubjects, 'mcorrect', mcorrect, 'condstat', 'off', 'cond2stat', 'off', ...
                                 'statistics', method, 'condnames', allconditions, 'cond2names', allgroups, 'clustname', STUDY.cluster(opt.clusters(index)).name, 'compnames', comp_names, ...
                                 'subject', opt.subject, 'valsunit', opt.unitx, 'vals', params.topotime, 'datatype', datatypestr, 'cond2group', params.plotgroups, 'condgroup', params.plotconditions);
                             
            handles = findall(0,'Type','Figure', 'Tag','tmp_curvetag');
            std_detachplots('','','data',erpdata2,'figtitles', {alltitlestmp{:}}','sbtitles',sbtitles,'handles', handles, 'filter',params.filter);
            axcopyflag = 0;

        end
%--------------------------------------------------------------------------
    end
    tmpgcf = gcf;
    set(tmpgcf,'name', ['Component ' datatypestr ] );
    if axcopyflag
          haxis = findall(tmpgcf,'type','axes');
        for i= 1: length(haxis)
            axcopy(haxis(i));
        end
    end 
end

% remove fields and ignore fields who are absent
% ----------------------------------------------
function s = myrmfield(s, f)

for index = 1:length(f)
    if isfield(s, f{index})
        s = rmfield(s, f{index});
    end
end

% convert to structure (but take into account cells)
% --------------------------------------------------
function s = mystruct(v)

for index=1:length(v)
    if iscell(v{index})
        v{index} = { v{index} };
    end
end
try
    s = struct(v{:});
catch, error('Parameter error'); end

% convert to structure (but take into account cells)
% --------------------------------------------------
function s = myfieldnames(v)

s = fieldnames(v);
if isfield(v, 'eeglab')
    s2 = fieldnames(v.eeglab);
    s = { s{:} s2{:} };
end
if isfield(v, 'fieldtrip')
    s3 = fieldnames(v.fieldtrip);
    for index=1:length(s3)
        s3{index} = [ 'fieldtrip' s3{index} ];
    end
    s = { s{:} s3{:} };
end
