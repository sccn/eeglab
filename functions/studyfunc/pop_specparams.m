% pop_specparams() - Set plotting and statistics parameters for computing
%                    STUDY component spectra.
% Usage:    
%       >> STUDY = pop_specparams(STUDY, 'key', 'val');   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% Statistics options:
%   'groupstats' - ['on'|'off'] Compute (or not) statistics across subject 
%                  groups {default: 'off'}
%   'condstats'  - ['on'|'off'] Compute (or not) statistics across data.
%                  conditions {default: 'off'}
%   'topofreq'   - [real] Plot Spectrum scalp maps at one specific freq. (Hz).
%                  A frequency range [min max] may also be defined (the 
%                  spectrum is then averaged over the interval) {default: []}
%   'statistics' - ['param'|'perm'|'bootstrap'] Type of statistics to compute
%                  'param' for parametric (t-test/anova); 'perm' for 
%                  permutation-based and 'bootstrap' for bootstrap 
%                  {default: 'param'}
%   'naccu'      - [integer] Number of surrogate averages to accumulate for
%                  surrogate statistics. For example, to test whether 
%                  use 'naccu' >= 200; for p < 0.001, use >= 2000. When 
%                  threshold (below) is not NaN and 'naccu' is too low, 
%                  'naccu' will be automatically updated (for now, from
%                  the command line only).
%   'threshold'  - [NaN|0.0x] Significance threshold. NaN will plot the 
%                  p-value themselves on a different figure. When possible, 
%                  significant latency regions are shown below the data.
%                  {default: NaN}
%   'mcorrect'   - ['fdr'|'none'] correction for multiple comparisons
%                  (threshold case only). 'fdr' uses false discovery rate.
%                  See the fdr function for more information. Defaut is
%                  'none'.
% Plot options:
%   'freqrange'  - [min max] spectral frequency range (in Hz) to plot. 
%                  {default: whole frequency range} .
%   'ylim'       - [mindB maxdB] spectral plotting limits in dB 
%                  {default: from data}
%   'plotgroups' - ['together'|'apart'] 'together' -> plot subject groups 
%                  on the same figure in different colors, else ('apart') on 
%                  different figures {default: 'apart'}
%   'plotconditions' - ['together'|'apart'] 'together' -> plot conditions 
%                  on the same figure in different colors, else ('apart') 
%                  on different figures. Note: keywords 'plotgroups' and 
%                  'plotconditions' cannot both be set to 'together'. 
%                  {default: 'apart'}
%
% See also: std_specplot()
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2006-

% Copyright (C) Arnaud Delorme, CERCO
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

% $Log: pop_specparams.m,v $
% Revision 1.20  2010/02/26 10:57:36  claire
% typo
%
% Revision 1.19  2010/02/24 10:52:36  arno
% Implemented new single trial statistics
%
% Revision 1.18  2009/11/11 00:28:53  arno
% New GUI format
%
% Revision 1.17  2009/05/31 04:37:38  arno
% Better GUI for Study plotting options
%
% Revision 1.16  2009/05/31 02:22:10  arno
% Adding FDR and bootstrap to all STUDY functions
%
% Revision 1.15  2007/10/25 18:50:07  nima
% ??s filled.
%
% Revision 1.14  2007/03/17 21:23:04  arno
% logical operator precedence
%
% Revision 1.13  2007/01/26 18:03:27  arno
% new topofreq option
%
% Revision 1.12  2006/11/23 00:59:36  arno
% removing psectral data if subject mean subtract
%
% Revision 1.11  2006/11/23 00:48:13  arno
% typo
%
% Revision 1.10  2006/11/23 00:47:41  arno
% subject mean buttont
%
% Revision 1.9  2006/11/10 01:29:49  arno
% GUI size
%
% Revision 1.8  2006/11/04 00:19:43  arno
% text
%
% Revision 1.7  2006/10/03 21:29:15  scott
% help msg edits.  Some ?? remain.  -sm
%
% Revision 1.6  2006/10/03 18:10:30  scott
% minor help edit
%
% Revision 1.5  2006/10/02 20:58:43  scott
% plotcond -> plotconditions
%
% Revision 1.4  2006/10/02 17:19:49  scott
% edited help msg for clarity and grammar. changed to English 'apart'
% changed 'plotgroup' to 'plotgroups'; in help, changed 'panel' to 'axis'
%
% Revision 1.3  2006/10/02 11:39:15  arno
% header documentation
%

function [ STUDY, com ] = pop_specparams(STUDY, varargin);

STUDY = default_params(STUDY);
STUDY.etc.specparams = pop_statparams(STUDY.etc.specparams, 'default');
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    enablegroup = fastif(length(STUDY.group)>1, 'on', 'off');
    enablecond  = fastif(length(STUDY.condition)>1, 'on', 'off');
    plotconditions    = fastif(strcmpi(STUDY.etc.specparams.plotconditions, 'together'), 1, 0);
    plotgroups   = fastif(strcmpi(STUDY.etc.specparams.plotgroups,'together'), 1, 0);
    submean      = fastif(strcmpi(STUDY.etc.specparams.subtractsubjectmean,'on'), 1, 0);
    vis = fastif(isnan(STUDY.etc.specparams.topofreq), 'off', 'on');
    
    uilist = { ...
        {'style' 'text'       'string' 'Frequency [low_Hz high_Hz]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.specparams.freqrange) 'tag' 'freqrange' } ...
        {'style' 'text'       'string' 'Plot limits [low high]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.specparams.ylim) 'tag' 'ylim' } ...
        {'style' 'text'       'string' 'Plot scalp map at freq. [Hz]' 'enable' vis } ...
        {'style' 'edit'       'string' num2str(STUDY.etc.specparams.topofreq) 'tag' 'topofreq' 'enable' vis } ...
        {} {} ...
        {} {'style' 'checkbox'   'string' 'Subtract individual subject mean spectrum' 'value' submean 'tag' 'submean' } ...
        {} {'style' 'checkbox'   'string' 'Plot conditions on the same panel' 'value' plotconditions 'enable' enablecond  'tag' 'plotconditions' } ...
        {} {'style' 'checkbox'   'string' 'Plot groups on the same panel' 'value' plotgroups 'enable' enablegroup 'tag' 'plotgroups' } };

    cbline = [0.07 1.1];
    otherline = [ 0.7 .5 0.6 .5];
    geometry = { otherline otherline cbline cbline cbline };
    
    [STUDY.etc.specparams res options] = pop_statparams(STUDY.etc.specparams, 'geometry' , geometry, 'uilist', uilist, ...
                                   'helpcom', 'pophelp(''std_specparams'')', 'enablegroup', enablegroup, ...
                                   'enablecond', enablecond, ...
                                   'title', 'Set spectrum plotting parameters -- pop_specparams()');

    if isempty(res), return; end;
    
    % decode inputs
    % -------------
    if res.plotgroups & res.plotconditions, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end;
    if res.submean   , res.submean    = 'on'; else res.submean    = 'off'; end;
    if res.plotgroups, res.plotgroups = 'together'; else res.plotgroups = 'apart'; end;
    if res.plotconditions , res.plotconditions  = 'together'; else res.plotconditions  = 'apart'; end;
    res.topofreq  = str2num( res.topofreq );
    res.freqrange = str2num( res.freqrange );
    res.ylim      = str2num( res.ylim );
    
    % build command call
    % ------------------
    if ~strcmpi( res.plotgroups, STUDY.etc.specparams.plotgroups), options = { options{:} 'plotgroups' res.plotgroups }; end;
    if ~strcmpi( res.plotconditions , STUDY.etc.specparams.plotconditions ), options = { options{:} 'plotconditions'  res.plotconditions  }; end;
    if ~strcmpi( res.submean   , STUDY.etc.specparams.subtractsubjectmean ), options = { options{:} 'subtractsubjectmean'  res.submean  }; end;
    if ~isequal(res.topofreq, STUDY.etc.specparams.topofreq),   options = { options{:} 'topofreq' res.topofreq }; end;
    if ~isequal(res.ylim, STUDY.etc.specparams.ylim),           options = { options{:} 'ylim' res.ylim      }; end;
    if ~isequal(res.freqrange, STUDY.etc.specparams.freqrange), options = { options{:} 'freqrange' res.freqrange }; end;
    if ~isempty(options)
        STUDY = pop_specparams(STUDY, options{:});
        com = sprintf('STUDY = pop_specparams(STUDY, %s);', vararg2str( options ));
    end;
else
    
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            STUDY.etc.specparams = setfield(STUDY.etc.specparams, varargin{index}, varargin{index+1});
        end;
    end;
end;

% scan clusters and channels to remove specdata info if freqrange has changed
% ----------------------------------------------------------
if ~isequal(STUDY.etc.specparams.freqrange, TMPSTUDY.etc.specparams.freqrange) | ...
    ~isequal(STUDY.etc.specparams.subtractsubjectmean, TMPSTUDY.etc.specparams.subtractsubjectmean)
    if isfield(STUDY.cluster, 'specdata')
        for index = 1:length(STUDY.cluster)
            STUDY.cluster(index).specdata  = [];
            STUDY.cluster(index).specfreqs = [];
        end;
    end;
    if isfield(STUDY.changrp, 'specdata')
        for index = 1:length(STUDY.changrp)
            STUDY.changrp(index).specdata  = [];
            STUDY.changrp(index).specfreqs = [];
        end;
    end;
end;

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'specparams'), STUDY.etc.specparams = []; end;
    if ~isfield(STUDY.etc.specparams, 'topofreq'),   STUDY.etc.specparams.topofreq = []; end;
    if ~isfield(STUDY.etc.specparams, 'freqrange'),  STUDY.etc.specparams.freqrange = []; end;
    if ~isfield(STUDY.etc.specparams, 'ylim'     ),  STUDY.etc.specparams.ylim      = []; end;
    if ~isfield(STUDY.etc.specparams, 'subtractsubjectmean' ), STUDY.etc.specparams.subtractsubjectmean  = 'off'; end;
    if ~isfield(STUDY.etc.specparams, 'plotgroups'), STUDY.etc.specparams.plotgroups = 'apart'; end;
    if ~isfield(STUDY.etc.specparams, 'plotconditions'),  STUDY.etc.specparams.plotconditions  = 'apart'; end;
