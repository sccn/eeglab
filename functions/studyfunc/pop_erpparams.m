% pop_erpparams() - Set plotting and statistics parameters for cluster ERP 
%                   plotting
% Usage:    
%   >> STUDY = pop_erpparams(STUDY, 'key', 'val');   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% Statistics options:
%   'groupstats  - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%   'condstats'  - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%   'topotime'   - [real] Plot ERP scalp maps at one specific latency (ms).
%                   A latency range [min max] may also be defined (the 
%                   ERP is then averaged over the interval) {default: []}
%   'statistics' - ['param'|'perm'|'bootstrap'] Type of statistics to compute
%                  'param' for parametric (t-test/anova); 'perm' for 
%                  permutation-based and 'bootstrap' for bootstrap 
%                  {default: 'param'}
%   'naccu'      - [integer] Number of surrogate averages to accumulate for
%                  surrogate statistics. For example, to test whether 
%                  p<0.01, use >=200. For p<0.001, use 'naccu' >=2000. 
%                  If a threshold (not NaN) is set below, and 'naccu' is 
%                  too low, it will be automatically reset. (This option 
%                  is now available only from the command line).
%   'threshold'  - [NaN|float<<1] Significance probability threshold. 
%                  NaN -> plot the p-values themselves on a different axis. 
%                  When possible, the significant time regions are indicated 
%                  below the data.
%   'mcorrect'   - ['fdr'|'none'] correction for multiple comparisons
%                  (threshold case only). 'fdr' uses false discovery rate.
%                  See the fdr function for more information. Defaut is
%                  'none'.
% Plot options:
%   'timerange'  - [min max] ERP plotting latency range in ms. 
%                  {default: the whole epoch}
%   'ylim'       - [min max] ERP limits in microvolts {default: from data}
%   'plotgroups' - ['together'|'apart'] 'together' -> plot subject groups 
%                  on the same axis in different colors, else ('apart') 
%                  on different axes. {default: 'apart'}
%   'plotconditions' - ['together'|'apart'] 'together' -> plot conditions 
%                  on the same axis in different colors, else ('apart') 
%                  on different axes. Note: Keywords 'plotgroups' and 
%                  'plotconditions' may not both be set to 'together'. 
%                  {default: 'together'}
% See also: std_erpplot()
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2006-

% Copyright (C) Arnaud Delorme, 2006
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
% Revision 1.21  2009/05/31 02:22:10  arno
% Adding FDR and bootstrap to all STUDY functions
%
% Revision 1.20  2007/11/21 16:48:29  arno
% remove ??
%
% Revision 1.19  2007/08/12 03:41:29  arno
% fixing topotime
%
% Revision 1.18  2007/08/12 03:38:15  arno
% change code for scalp latency plotting
%
% Revision 1.17  2007/08/12 03:36:20  arno
% typo
%
% Revision 1.16  2007/03/20 03:35:14  arno
% fixed topotime
%
% Revision 1.15  2007/03/20 03:08:59  arno
% fix filter history option
%
% Revision 1.14  2007/03/17 21:20:06  arno
% logical operator precedence
%
% Revision 1.13  2007/01/26 18:02:16  arno
% new topotime option
%
% Revision 1.12  2006/11/22 20:04:40  arno
% same
%
% Revision 1.11  2006/11/22 20:03:35  arno
% filter filed
% field
%
% Revision 1.10  2006/11/10 01:30:21  arno
% gui size
%
% Revision 1.9  2006/11/10 01:28:28  arno
% GUI size
%
% Revision 1.8  2006/11/04 00:20:25  arno
% edit text
%
% Revision 1.7  2006/10/03 21:20:25  scott
% helkp msg edits - see remaining few ??  -sm
%
% Revision 1.6  2006/10/03 18:09:34  scott
% minor help edit
%
% Revision 1.5  2006/10/02 21:57:51  scott
% plotcond -> plotconditions
%
% Revision 1.4  2006/10/02 17:09:09  scott
% edited help message for clarity and grammar. NOTE: changed argument
% 'appart'  to English 'apart'. ALSO changed keyword 'plotgroup' to
% more grammatical 'plotgroups' throughout. (Note: did not change
% 'plotcond' to 'plotconds' since this abbreviation feels awkward).
% Similar changes may be needed in other functions. -sm
%
% Revision 1.3  2006/10/02 11:38:06  arno
% header documentation
%

function [ STUDY, com ] = pop_erpparams(STUDY, varargin);

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    enablegroup = fastif(length(STUDY.group)>1, 'on', 'off');
    enablecond  = fastif(length(STUDY.condition)>1, 'on', 'off');
    threshstr   = fastif(isnan(STUDY.etc.erpparams.threshold),'', num2str(STUDY.etc.erpparams.threshold));
    plotconditions    = fastif(strcmpi(STUDY.etc.erpparams.plotconditions, 'together'), 1, 0);
    plotgroups  = fastif(strcmpi(STUDY.etc.erpparams.plotgroups,'together'), 1, 0);
    condstats    = fastif(strcmpi(STUDY.etc.erpparams.condstats, 'on'), 1, 0);
    groupstats   = fastif(strcmpi(STUDY.etc.erpparams.groupstats,'on'), 1, 0);
    mcorrect     = fastif(strcmpi(STUDY.etc.erpparams.mcorrect,  'fdr'), 1, 0);
    if strcmpi(STUDY.etc.erpparams.statistics,'param'),    statval = 1;
    elseif strcmpi(STUDY.etc.erpparams.statistics,'perm'), statval = 2;
    else                                                   statval = 3;
    end;
    vis = fastif(isnan(STUDY.etc.erpparams.topotime), 'off', 'on');
    
    uilist = { ...
        {'style' 'text'       'string' 'Time range in ms [low high]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.timerange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Plot limits in uV [low high]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.ylim) 'tag' 'ylim' } ...
        {'style' 'text'       'string' 'Plot scalp map at latency [ms]' 'enable' vis } ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.topotime) 'tag' 'topotime' 'enable' vis } ...
        {'style' 'text'       'string' 'Display filter in Hz [high]' } ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.filter) 'tag' 'filter' } ...
        {} {'style' 'checkbox'   'string' '' 'value' plotconditions 'enable' enablecond  'tag' 'plotconditions' } ...
        {'style' 'text'       'string' 'Plot conditions on the same panel' 'enable' enablecond } ...
        {} {'style' 'checkbox'   'string' '' 'value' plotgroups 'enable' enablegroup 'tag' 'plotgroups' } ...
        {'style' 'text'       'string' 'Plot groups on the same panel' 'enable' enablegroup } ...
        {} ...
        {'style' 'text'       'string' 'Statistics'} ...
        {'style' 'popupmenu'  'string' 'Parametric|Permutations|Bootstrap' 'tag' 'statistics' 'value' statval 'listboxtop' statval } ...
        {'style' 'text'       'string' 'Threshold (p<)' } ...
        {'style' 'edit'       'string' threshstr 'tag' 'threshold' } ...
        {} {'style' 'checkbox'   'string' '' 'value' condstats  'enable' enablecond  'tag' 'condstats' } ...
        {'style' 'text'       'string' 'Compute condition statistics' 'enable' enablecond} ...
        {} {'style' 'checkbox'   'string' '' 'value' groupstats 'enable' enablegroup 'tag' 'groupstats' } ...
        {'style' 'text'       'string' 'Compute group statistics' 'enable' enablegroup } ...
        {} {'style' 'checkbox'  'string' '' 'value' mcorrect 'tag' 'mcorrect' } ... 
        {'style' 'text' 'string' 'Use False Discovery Rate to correct for multiple comparisons' } };
    
    geometry = { [ 1 .5 1 .5] [1 0.5  1 0.5] [0.1 0.1 1] [0.1 0.1 1] [1] [.5 1 1 .5] [0.1 0.1 1] [0.1 0.1 1] [0.1 0.1 1] };
    
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, ...
                                   'helpcom', 'pophelp(''std_erpparams'')', ...
                                   'title', 'Set parameters for plotting ERPs -- pop_erpparams()');

    if isempty(res), return; end;
    
    % decode inputs
    % -------------
    if res.plotgroups & res.plotconditions, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end;
    if res.groupstats, res.groupstats = 'on'; else res.groupstats = 'off'; end;
    if res.condstats , res.condstats  = 'on'; else res.condstats  = 'off'; end;
    if res.mcorrect,   res.mcorrect   = 'fdr'; else res.mcorrect  = 'none'; end;
    if res.plotgroups, res.plotgroups = 'together'; else res.plotgroups = 'apart'; end;
    if res.plotconditions , res.plotconditions  = 'together'; else res.plotconditions  = 'apart'; end;
    res.topotime  = str2num( res.topotime );
    res.timerange = str2num( res.timerange );
    res.ylim      = str2num( res.ylim );
    res.threshold = str2num( res.threshold );
    res.filter    = str2num( res.filter );
    if isempty(res.threshold),res.threshold = NaN; end;
    if res.statistics == 1, res.statistics  = 'param'; 
    elseif res.statistics == 2, res.statistics  = 'perm'; 
    else res.statistics  = 'bootstrap'; 
    end;
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( char(res.filter), char(STUDY.etc.erpparams.filter)), options = { options{:} 'filter' res.filter }; end;
    if ~strcmpi( res.plotgroups, STUDY.etc.erpparams.plotgroups), options = { options{:} 'plotgroups' res.plotgroups }; end;
    if ~strcmpi( res.plotconditions , STUDY.etc.erpparams.plotconditions ), options = { options{:} 'plotconditions'  res.plotconditions  }; end;
    if ~strcmpi( res.groupstats, STUDY.etc.erpparams.groupstats), options = { options{:} 'groupstats' res.groupstats }; end;
    if ~strcmpi( res.condstats , STUDY.etc.erpparams.condstats ), options = { options{:} 'condstats'  res.condstats  }; end;
    if ~strcmpi( res.mcorrect,   STUDY.etc.erpparams.mcorrect),   options = { options{:} 'mcorrect' res.mcorrect }; end;
    if ~strcmpi( res.statistics, STUDY.etc.erpparams.statistics ), options = { options{:} 'statistics' res.statistics }; end;
    if ~isequal(res.ylim     , STUDY.etc.erpparams.ylim),      options = { options{:} 'ylim' res.ylim       }; end;
    if ~isequal(res.timerange, STUDY.etc.erpparams.timerange), options = { options{:} 'timerange' res.timerange }; end;
    if (all(isnan(res.topotime)) & all(~isnan(STUDY.etc.erpparams.topotime))) | ...
            (all(~isnan(res.topotime)) & all(isnan(STUDY.etc.erpparams.topotime))) | ...
                (all(~isnan(res.topotime)) & ~isequal(res.topotime, STUDY.etc.erpparams.topotime))
        options = { options{:} 'topotime' res.topotime }; 
    end;
    if (isnan(res.threshold) & ~isnan(STUDY.etc.erpparams.threshold)) | ...
            (~isnan(res.threshold) & isnan(STUDY.etc.erpparams.threshold)) | ...
                (~isnan(res.threshold) & res.threshold ~= STUDY.etc.erpparams.threshold)
                options = { options{:} 'threshold' res.threshold }; 
    end;
    if ~isempty(options)
        STUDY = pop_erpparams(STUDY, options{:});
        com = sprintf('STUDY = pop_erpparams(STUDY, %s);', vararg2str( options ));
    end;
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            STUDY.etc.erpparams = setfield(STUDY.etc.erpparams, varargin{index}, varargin{index+1});
        end;
    end;
end;

% scan clusters and channels to remove erpdata info if timerange has changed
% ----------------------------------------------------------
if ~isequal(STUDY.etc.erpparams.timerange, TMPSTUDY.etc.erpparams.timerange)
    if isfield(STUDY.cluster, 'erpdata')
        for index = 1:length(STUDY.cluster)
            STUDY.cluster(index).erpdata  = [];
            STUDY.cluster(index).erptimes = [];
        end;
    end;
    if isfield(STUDY.changrp, 'erpdata')
        for index = 1:length(STUDY.changrp)
            STUDY.changrp(index).erpdata  = [];
            STUDY.changrp(index).erptimes = [];
        end;
    end;
end;

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'erpparams'), STUDY.etc.erpparams = []; end;
    if ~isfield(STUDY.etc.erpparams, 'topotime'),    STUDY.etc.erpparams.topotime = []; end;
    if ~isfield(STUDY.etc.erpparams, 'filter'),     STUDY.etc.erpparams.filter = []; end;
    if ~isfield(STUDY.etc.erpparams, 'timerange'),  STUDY.etc.erpparams.timerange = []; end;
    if ~isfield(STUDY.etc.erpparams, 'ylim'     ),  STUDY.etc.erpparams.ylim      = []; end;
    if ~isfield(STUDY.etc.erpparams, 'statistics'), STUDY.etc.erpparams.statistics = 'param'; end;
    if ~isfield(STUDY.etc.erpparams, 'groupstats'),  STUDY.etc.erpparams.groupstats = 'off'; end;
    if ~isfield(STUDY.etc.erpparams, 'condstats' ),  STUDY.etc.erpparams.condstats  = 'off'; end;
    if ~isfield(STUDY.etc.erpparams, 'mcorrect' ),   STUDY.etc.erpparams.mcorrect   = 'none'; end;
    if ~isfield(STUDY.etc.erpparams, 'threshold' ),  STUDY.etc.erpparams.threshold = NaN; end;
    if ~isfield(STUDY.etc.erpparams, 'plotgroups') , STUDY.etc.erpparams.plotgroups = 'apart'; end;
    if ~isfield(STUDY.etc.erpparams, 'naccu') ,      STUDY.etc.erpparams.naccu     = []; end;
    if ~isfield(STUDY.etc.erpparams, 'plotconditions') ,  STUDY.etc.erpparams.plotconditions  = 'apart'; end;

