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

% $Log: pop_erpparams.m,v $
% Revision 1.24  2010/02/24 10:52:36  arno
% Implemented new single trial statistics
%
% Revision 1.23  2009/11/11 00:28:53  arno
% New GUI format
%
% Revision 1.22  2009/05/31 04:37:38  arno
% Better GUI for Study plotting options
%
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
STUDY.etc.erpparams = pop_statparams(STUDY.etc.erpparams, 'default');
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    enablecond  = fastif(length(STUDY.design(STUDY.currentdesign).indvar1)>1, 'on', 'off');
    enablegroup = fastif(length(STUDY.design(STUDY.currentdesign).indvar2)>1, 'on', 'off');
    plotconditions    = fastif(strcmpi(STUDY.etc.erpparams.plotconditions, 'together'), 1, 0);
    plotgroups  = fastif(strcmpi(STUDY.etc.erpparams.plotgroups,'together'), 1, 0);
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
        {} {'style' 'checkbox'   'string' 'Plot conditions on the same panel' 'value' plotconditions 'enable' enablecond  'tag' 'plotconditions' } ...
        {} {'style' 'checkbox'   'string' 'Plot groups on the same panel' 'value' plotgroups 'enable' enablegroup 'tag' 'plotgroups' } };
    
    cbline = [0.07 1.1];
    otherline = [ 0.7 .5 0.6 .5];
    geometry = { otherline otherline cbline cbline };
    
    [STUDY.etc.erpparams res options] = pop_statparams(STUDY.etc.erpparams, 'geometry' , geometry, 'uilist', uilist, ...
                                   'helpcom', 'pophelp(''std_erpparams'')', 'enablegroup', enablegroup, ...
                                   'enablecond', enablecond, ...
                                   'title', 'Set ERP plotting parameters -- pop_erpparams()');
                               
    if isempty(res), return; end;
    
    % decode inputs
    % -------------
    if res.plotgroups & res.plotconditions, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end;
    if res.plotgroups, res.plotgroups = 'together'; else res.plotgroups = 'apart'; end;
    if res.plotconditions , res.plotconditions  = 'together'; else res.plotconditions  = 'apart'; end;
    res.topotime  = str2num( res.topotime );
    res.timerange = str2num( res.timerange );
    res.ylim      = str2num( res.ylim );
    res.filter    = str2num( res.filter );
    
    % build command call
    % ------------------
    if ~strcmpi( char(res.filter), char(STUDY.etc.erpparams.filter)), options = { options{:} 'filter' res.filter }; end;
    if ~strcmpi( res.plotgroups, STUDY.etc.erpparams.plotgroups), options = { options{:} 'plotgroups' res.plotgroups }; end;
    if ~strcmpi( res.plotconditions , STUDY.etc.erpparams.plotconditions ), options = { options{:} 'plotconditions'  res.plotconditions  }; end;
    if ~isequal(res.ylim     , STUDY.etc.erpparams.ylim),      options = { options{:} 'ylim' res.ylim       }; end;
    if ~isequal(res.timerange, STUDY.etc.erpparams.timerange), options = { options{:} 'timerange' res.timerange }; end;
    if (all(isnan(res.topotime)) & all(~isnan(STUDY.etc.erpparams.topotime))) | ...
            (all(~isnan(res.topotime)) & all(isnan(STUDY.etc.erpparams.topotime))) | ...
                (all(~isnan(res.topotime)) & ~isequal(res.topotime, STUDY.etc.erpparams.topotime))
        options = { options{:} 'topotime' res.topotime }; 
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
    if ~isfield(STUDY.etc.erpparams, 'plotgroups') , STUDY.etc.erpparams.plotgroups = 'apart'; end;
    if ~isfield(STUDY.etc.erpparams, 'plotconditions') ,  STUDY.etc.erpparams.plotconditions  = 'apart'; end;

