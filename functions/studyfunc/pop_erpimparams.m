% pop_erpimparams() - Set plotting and statistics parameters for 
%                    computing and plotting STUDY mean ERPimages and measure 
%                    statistics. Settings are stored within the STUDY 
%                    structure (STUDY.etc.erpimparams) which is used
%                    whenever plotting is performed by the function
%                    std_erpimage().
% Usage:    
%   >> STUDY = pop_erpimparams(STUDY, 'key', 'val', ...);   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% Statistics options:
%  'groupstats'   - ['on'|'off'] Compute statistics across subject 
%                  groups {default: 'off'}
%  'condstats'    - ['on'|'off'] Compute statistics across data 
%                  conditions {default: 'off'}
%  'statistics'  - ['param'|'perm'|'bootstrap'] Type of statistics to compute
%                  'param' for parametric (t-test/anova); 'perm' for 
%                  permutation-based and 'bootstrap' for bootstrap 
%                  {default: 'param'}
%  'statmode'    - ['subjects'|'trials'] 'subjects' {default}
%                  -> statistics are computed across condition mean 
%                  erpims|ITCs of the single subjects. 
%                  'trials' -> single-trial 'erpim' transforms 
%                  for all subjects are pooled.  This requires that 
%                  they were saved to disk using std_erpim() option 
%                  'savetrials', 'on'. Note, however, that the 
%                  single-trial erpims may occupy several GB of disk 
%                  space, and that computation of statistics may 
%                  require a large amount of RAM.
%  'naccu'       - [integer] Number of surrogate data averages to use in
%                  surrogate statistics. For instance, if p<0.01, 
%                  use naccu>200. For p<0.001, naccu>2000. If a 'threshold'
%                  (not NaN) is set below and 'naccu' is too low, it will
%                  be automatically increased. (This keyword is currently
%                  only modifiable from the command line, not from the gui). 
%  'threshold'   - [NaN|alpha] Significance threshold (0<alpha<<1). Value 
%                  NaN will plot p-values for each time and/or frequency
%                  on a different axis. If alpha is used, significant time
%                  and/or frequency regions will be indicated either on
%                  a separate axis or (whenever possible) along with the
%                  data {default: NaN}
%   'mcorrect'   - ['fdr'|'none'] correction for multiple comparisons
%                  (threshold case only). 'fdr' uses false discovery rate.
%                  See the fdr function for more information. Defaut is
%                  'none'.
%  'subbaseline' - ['on'|'off'] subtract the same baseline across conditions 
%                  for erpim (not ITC). When datasets with different conditions
%                  are recorded simultaneously, a common baseline spectrum 
%                  should be used. Note that this also affects the 
%                  results of statistics {default: 'on'}
%  'maskdata'    - ['on'|'off'] when threshold is not NaN, and 'groupstats'
%                  or 'condstats' (above) are 'off', masks the data 
%                  for significance.
%
% erpim/ITC image plotting options:
%  'timerange'   - [min max] erpim/ITC plotting latency range in ms. 
%                  {default: the whole output latency range}.
%  'trialrange'  - [min max] erpim/ITC plotting frequency range in ms. 
%                  {default: the whole output frequency range}
%  'topotime'    - [float] plot scalp map at specific time. A time range may
%                  also be provide and the erpim will be averaged over the
%                  given time range. Requires 'topofreq' below to be set.
%  'topotrial'  - [float] plot scalp map at specific trial in ERPimage. As 
%                  above a trial range may also be provided.
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

function [ STUDY, com ] = pop_erpimparams(STUDY, varargin);

STUDY = default_params(STUDY);
STUDY.etc.erpimparams = pop_statparams(STUDY.etc.erpimparams, 'default');
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    vis = fastif(isnan(STUDY.etc.erpimparams.topotime), 'off', 'on');
    uilist = { ...
        {'style' 'text'       'string' 'Time range in ms [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.timerange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Plot scalp map at time [ms]' 'visible' vis} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.topotime) 'tag' 'topotime' 'visible' vis } ...
        {'style' 'text'       'string' 'Trial range [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.trialrange) 'tag' 'trialrange' } ...
        {'style' 'text'       'string' 'Plot scalp map at trial(s)' 'visible' vis} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.topotrial) 'tag' 'topotrial' 'visible' vis } ...
        {'style' 'text'       'string' 'Color limits [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.colorlimits) 'tag' 'colorlimits' } ...
        {'style' 'text'       'string' '' } ...
        {'style' 'text'       'string' '' } };
    
    cbline = [0.07 1.1];
    otherline = [ 0.7 .5 0.6 .5];
    geometry = { otherline otherline otherline };
    enablecond  = fastif(length(STUDY.design(STUDY.currentdesign).variable(1).value)>1, 'on', 'off');
    enablegroup = fastif(length(STUDY.design(STUDY.currentdesign).variable(2).value)>1, 'on', 'off');
    
    [STUDY.etc.erpimparams res options] = pop_statparams(STUDY.etc.erpimparams, 'geometry' , geometry, 'uilist', uilist, ...
                                   'helpcom', 'pophelp(''std_erpimparams'')', 'enablegroup', enablegroup, ...
                                   'enablecond', enablecond, 'enablesingletrials', 'off', ...
                                   'title', 'Set ERP-image plotting parameters -- pop_erpimparams()');
                               
    if isempty(res), return; end;
    
    % decode input
    % ------------
    res.topotime    = str2num( res.topotime );
    res.topotrial   = str2num( res.topotrial );
    res.timerange   = str2num( res.timerange );
    res.freqrange   = str2num( res.trialrange );
    res.colorlimits = str2num( res.colorlimits );
    
    % build command call
    % ------------------
    if ~isequal(res.topotime  ,  STUDY.etc.erpimparams.topotime),    options = { options{:} 'topotime'    res.topotime    }; end;
    if ~isequal(res.topotrial,   STUDY.etc.erpimparams.topotrial),   options = { options{:} 'topotrial'   res.topotrial   }; end;
    if ~isequal(res.timerange ,  STUDY.etc.erpimparams.timerange),   options = { options{:} 'timerange'   res.timerange   }; end;
    if ~isequal(res.trialrange,  STUDY.etc.erpimparams.trialrange),  options = { options{:} 'trialrange'  res.trialrange  }; end;
    if ~isequal(res.colorlimits, STUDY.etc.erpimparams.colorlimits), options = { options{:} 'colorlimits' res.colorlimits }; end;
    if ~isempty(options)
        STUDY = pop_erpimparams(STUDY, options{:});
        com = sprintf('STUDY = pop_erpimparams(STUDY, %s);', vararg2str( options ));
    end;
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            STUDY.etc.erpimparams = setfield(STUDY.etc.erpimparams, varargin{index}, varargin{index+1});
        end;
    end;
end;

% scan clusters and channels to remove erpimdata info if timerange etc. have changed
% ---------------------------------------------------------------------------------
if ~isequal(STUDY.etc.erpimparams.timerange, TMPSTUDY.etc.erpimparams.timerange) | ... 
    ~isequal(STUDY.etc.erpimparams.trialrange, TMPSTUDY.etc.erpimparams.trialrange)
    if isfield(STUDY.cluster, 'erpimdata')
        for index = 1:length(STUDY.cluster)
            STUDY.cluster(index).erpimdata   = [];
            STUDY.cluster(index).erpimtimes  = [];
            STUDY.cluster(index).erpimtrials = [];
            STUDY.cluster(index).erpimevents = [];
        end;
    end;
    if isfield(STUDY.changrp, 'erpimdata')
        for index = 1:length(STUDY.changrp)
            STUDY.changrp(index).erpimdata   = [];
            STUDY.changrp(index).erpimtimes  = [];
            STUDY.changrp(index).erpimtrials = [];
            STUDY.changrp(index).erpimevents = [];
        end;
    end;
end;

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'erpimparams'), STUDY.etc.erpimparams = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'erpimageopt'),  STUDY.etc.erpimparams.erpimageopt = {}; end;
    if ~isfield(STUDY.etc.erpimparams, 'sorttype'   ),  STUDY.etc.erpimparams.sorttype    = ''; end;
    if ~isfield(STUDY.etc.erpimparams, 'sortwin'    ),  STUDY.etc.erpimparams.sortwin     = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'sortfield'  ),  STUDY.etc.erpimparams.sortfield   = 'latency'; end;
    if ~isfield(STUDY.etc.erpimparams, 'statistics' ),  STUDY.etc.erpimparams.statistics  = 'param'; end;
    if ~isfield(STUDY.etc.erpimparams, 'groupstats' ),  STUDY.etc.erpimparams.groupstats  = 'off'; end;
    if ~isfield(STUDY.etc.erpimparams, 'condstats'  ),  STUDY.etc.erpimparams.condstats   = 'off'; end;
    if ~isfield(STUDY.etc.erpimparams, 'statmode'   ),  STUDY.etc.erpimparams.statmode    = 'param'; end;
    if ~isfield(STUDY.etc.erpimparams, 'threshold'  ),  STUDY.etc.erpimparams.threshold   = NaN; end;
    if ~isfield(STUDY.etc.erpimparams, 'mcorrect') ,    STUDY.etc.erpimparams.mcorrect    = 'none'; end;
    if ~isfield(STUDY.etc.erpimparams, 'naccu') ,       STUDY.etc.erpimparams.naccu       = []; end;

    if ~isfield(STUDY.etc.erpimparams, 'rmcomps'    ),  STUDY.etc.erpimparams.rmcomps     = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'interp'     ),  STUDY.etc.erpimparams.interp      = []; end;
    
    if ~isfield(STUDY.etc.erpimparams, 'timerange'  ),  STUDY.etc.erpimparams.timerange   = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'trialrange' ),  STUDY.etc.erpimparams.trialrange  = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'topotime'   ),  STUDY.etc.erpimparams.topotime    = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'topotrial' ),   STUDY.etc.erpimparams.topotrial   = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'colorlimits'),  STUDY.etc.erpimparams.colorlimits = []; end;

    if ~isfield(STUDY.etc.erpimparams, 'concatenate'),  STUDY.etc.erpimparams.concatenate = 'off'; end;
    if ~isfield(STUDY.etc.erpimparams, 'nlines'),       STUDY.etc.erpimparams.nlines      = 20; end;
    if ~isfield(STUDY.etc.erpimparams, 'smoothing'),    STUDY.etc.erpimparams.smoothing   = 10; end;
