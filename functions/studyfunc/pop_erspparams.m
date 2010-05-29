% pop_erspparams() - Set plotting and statistics parameters for 
%                    computing and plotting STUDY mean (and optionally 
%                    single-trial) ERSP and ITC measures and measure 
%                    statistics. Settings are stored within the STUDY 
%                    structure (STUDY.etc.erspparams) which is used
%                    whenever plotting is performed by the function
%                    std_specplot().
% Usage:    
%   >> STUDY = pop_erspparams(STUDY, 'key', 'val', ...);   
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
%                  ERSPs|ITCs of the single subjects. 
%                  'trials' -> single-trial 'ERSP' transforms 
%                  for all subjects are pooled.  This requires that 
%                  they were saved to disk using std_ersp() option 
%                  'savetrials', 'on'. Note, however, that the 
%                  single-trial ERSPs may occupy several GB of disk 
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
%                  for ERSP (not ITC). When datasets with different conditions
%                  are recorded simultaneously, a common baseline spectrum 
%                  should be used. Note that this also affects the 
%                  results of statistics {default: 'on'}
%  'maskdata'    - ['on'|'off'] when threshold is not NaN, and 'groupstats'
%                  or 'condstats' (above) are 'off', masks the data 
%                  for significance.
%
% ERSP/ITC image plotting options:
%  'timerange'   - [min max] ERSP/ITC plotting latency range in ms. 
%                  {default: the whole output latency range}.
%  'freqrange'   - [min max] ERSP/ITC plotting frequency range in ms. 
%                  {default: the whole output frequency range}
%  'ersplim'     - [mindB maxdB] ERSP plotting limits in dB 
%                  {default: from [ERSPmin,ERSPmax]}
%  'itclim'      - [minitc maxitc] ITC plotting limits (range: [0,1]) 
%                  {default: from [0,ITC data max]}
%  'topotime'    - [float] plot scalp map at specific time. A time range may
%                  also be provide and the ERSP will be averaged over the
%                  given time range. Requires 'topofreq' below to be set.
%  'topofreqs'   - [float] plot scalp map at specific frequencies. As above
%                  a frequency range may also be provided.
%
% See also: std_erspplot(), std_itcplot()
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

function [ STUDY, com ] = pop_erspparams(STUDY, varargin);

STUDY = default_params(STUDY);
STUDY.etc.erspparams = pop_statparams(STUDY.etc.erspparams, 'default');
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    subbaseline = fastif(strcmpi(STUDY.etc.erspparams.subbaseline,'on'), 1, 0);
    vis = fastif(isnan(STUDY.etc.erspparams.topotime), 'off', 'on');
    
    uilist = { ...
        {'style' 'text'       'string' 'Time range in ms [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erspparams.timerange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Plot scalp map at time [ms]' 'visible' vis} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erspparams.topotime) 'tag' 'topotime' 'visible' vis } ...
        {'style' 'text'       'string' 'Freq. range in Hz [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erspparams.freqrange) 'tag' 'freqrange' } ...
        {'style' 'text'       'string' 'Plot scalp map at freq. [Hz]' 'visible' vis} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erspparams.topofreq) 'tag' 'topofreq' 'visible' vis } ...
        {'style' 'text'       'string' 'Power limits in dB [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erspparams.ersplim) 'tag' 'ersplim' } ...
        {'style' 'text'       'string' 'ITC limit (0-1) [High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erspparams.itclim) 'tag' 'itclim' } ...
        {} {'style' 'checkbox'   'string' 'Compute ERSP baseline across conditions' 'value' subbaseline 'tag' 'subbaseline' }  };
    
    cbline = [0.07 1.1];
    otherline = [ 0.7 .5 0.6 .5];
    geometry = { otherline otherline otherline cbline };
    enablecond  = fastif(length(STUDY.design(STUDY.currentdesign).variable(1).value)>1, 'on', 'off');
    enablegroup = fastif(length(STUDY.design(STUDY.currentdesign).variable(2).value)>1, 'on', 'off');
    
    [STUDY.etc.erspparams res options] = pop_statparams(STUDY.etc.erspparams, 'geometry' , geometry, 'uilist', uilist, ...
                                   'helpcom', 'pophelp(''std_erspparams'')', 'enablegroup', enablegroup, ...
                                   'enablecond', enablecond, ...
                                   'title', 'Set ERSP|ITC plotting parameters -- pop_erspparams()');
                               
    if isempty(res), return; end;
    
    % decode input
    % ------------
    if res.subbaseline, res.subbaseline = 'on'; else res.subbaseline = 'off'; end;
    res.topotime  = str2num( res.topotime );
    res.topofreq  = str2num( res.topofreq );
    res.timerange = str2num( res.timerange );
    res.freqrange = str2num( res.freqrange );
    res.ersplim   = str2num( res.ersplim );
    res.itclim    = str2num( res.itclim );
    
    % build command call
    % ------------------
    if ~strcmpi( res.subbaseline , STUDY.etc.erspparams.subbaseline ), options = { options{:} 'subbaseline' res.subbaseline }; end;
    if ~isequal(res.topotime , STUDY.etc.erspparams.topotime),  options = { options{:} 'topotime'   res.topotime  }; end;
    if ~isequal(res.topofreq , STUDY.etc.erspparams.topofreq),  options = { options{:} 'topofreq'   res.topofreq  }; end;
    if ~isequal(res.ersplim  , STUDY.etc.erspparams.ersplim),   options = { options{:} 'ersplim'   res.ersplim   }; end;
    if ~isequal(res.itclim   , STUDY.etc.erspparams.itclim),    options = { options{:} 'itclim'    res.itclim    }; end;
    if ~isequal(res.timerange, STUDY.etc.erspparams.timerange), options = { options{:} 'timerange' res.timerange }; end;
    if ~isequal(res.freqrange, STUDY.etc.erspparams.freqrange), options = { options{:} 'freqrange' res.freqrange }; end;
    if ~isempty(options)
        STUDY = pop_erspparams(STUDY, options{:});
        com = sprintf('STUDY = pop_erspparams(STUDY, %s);', vararg2str( options ));
    end;
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            STUDY.etc.erspparams = setfield(STUDY.etc.erspparams, varargin{index}, varargin{index+1});
        end;
    end;
end;

% scan clusters and channels to remove erspdata info if timerange etc. have changed
% ---------------------------------------------------------------------------------
if ~isequal(STUDY.etc.erspparams.timerange, TMPSTUDY.etc.erspparams.timerange) | ... 
    ~isequal(STUDY.etc.erspparams.freqrange, TMPSTUDY.etc.erspparams.freqrange) | ... 
    ~isequal(STUDY.etc.erspparams.subbaseline, TMPSTUDY.etc.erspparams.subbaseline)
    if isfield(STUDY.cluster, 'erspdata')
        for index = 1:length(STUDY.cluster)
            STUDY.cluster(index).topotime  = [];
            STUDY.cluster(index).topofreq  = [];
            STUDY.cluster(index).erspdata  = [];
            STUDY.cluster(index).erspbase  = [];
            STUDY.cluster(index).ersptimes = [];
            STUDY.cluster(index).erspfreqs = [];
            STUDY.cluster(index).itcdata  = [];
            STUDY.cluster(index).itctimes = [];
            STUDY.cluster(index).itcfreqs = [];
        end;
    end;
    if isfield(STUDY.changrp, 'erspdata')
        for index = 1:length(STUDY.changrp)
            STUDY.changrp(index).topotime  = [];
            STUDY.changrp(index).topofreq  = [];
            STUDY.changrp(index).erspdata  = [];
            STUDY.changrp(index).erspbase  = [];
            STUDY.changrp(index).ersptimes = [];
            STUDY.changrp(index).erspfreqs = [];
            STUDY.changrp(index).itcdata  = [];
            STUDY.changrp(index).itctimes = [];
            STUDY.changrp(index).itcfreqs = [];
        end;
    end;
end;

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'erspparams'), STUDY.etc.erspparams = []; end;
    if ~isfield(STUDY.etc.erspparams, 'topotime'),     STUDY.etc.erspparams.topotime = []; end;
    if ~isfield(STUDY.etc.erspparams, 'topofreq'),     STUDY.etc.erspparams.topofreq = []; end;
    if ~isfield(STUDY.etc.erspparams, 'timerange'),    STUDY.etc.erspparams.timerange = []; end;
    if ~isfield(STUDY.etc.erspparams, 'freqrange'),    STUDY.etc.erspparams.freqrange = []; end;
    if ~isfield(STUDY.etc.erspparams, 'ersplim' ),     STUDY.etc.erspparams.ersplim   = []; end;
    if ~isfield(STUDY.etc.erspparams, 'itclim' ),      STUDY.etc.erspparams.itclim    = []; end;
    if ~isfield(STUDY.etc.erspparams, 'subbaseline' ),  STUDY.etc.erspparams.subbaseline = 'on';  end;

