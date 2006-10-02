% pop_specparams() - Set plotting and statistics parameters for STUDY spectra.
%
% Usage:    
%       >> STUDY = pop_specparams(STUDY, 'key', 'val');   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% Statistics options:
%   'statgroup'  - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%   'statcond'   - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%   'statistics' - ['param'|'perm'] Type of statistics to use 'param' for
%                  parametric and 'perm' for permutation-based statistics. 
%                  {default: 'param'}
%   'naccu'      - [integer] Number of surroage averages to accumulate for
%                  permutation-based statistics. For instance for p < 0.01,
%                  use 'naccu' >= 200; for p < 0.001, use >= 2000. When 
%                  threshold (below) is not NaN and 'naccu' is too low, 
%                  'naccu' will be automatically updated (for now, from
%                  the command line only).
%   'threshold'  - [NaN|0.0x] Significance threshold. NaN will plot the 
%                  p-value itself on a different axis. When possible, the
%                  significant latency regions are shown below the data.
% Plot options:
%   'freqrange'  - [min max] spectral frequency range (in Hz) to plot. 
%                  {default: whole frequency range} .
%   'ylim'       - [min max] spectral plotting limits {default: from data}
%   'plotgroups'  - ['together'|'apart'] 'together' -> plot subject groups 
%                  on the same axis in different colors, else ('apart') on 
%                  different panels.
%   'plotconditions' - ['together'|'apart'] 'together' -> plot conditions 
%                  on the same axis in different colors, else ('apart') 
%                  on different panels. Note: Keywords 'plotgroups' and 
%                  'plotconditions' cannot both be set to 'together'. 
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

% $Log: not supported by cvs2svn $
% Revision 1.4  2006/10/02 17:19:49  scott
% edited help msg for clarity and grammar. changed 'appart'  to English 'apart'
% changed 'plotgroup' to 'plotgroups'; in help, changed 'panel' to 'axis'
%
% Revision 1.3  2006/10/02 11:39:15  arno
% header documentation
%

function [ STUDY, com ] = pop_specparams(STUDY, varargin);

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    enablegroup = fastif(length(STUDY.group)>1, 'on', 'off');
    enablecond  = fastif(length(STUDY.condition)>1, 'on', 'off');
    threshstr   = fastif(isnan(STUDY.etc.specparams.threshold),'', num2str(STUDY.etc.specparams.threshold));
    plotconditions    = fastif(strcmpi(STUDY.etc.specparams.plotconditions, 'together'), 1, 0);
    plotgroups   = fastif(strcmpi(STUDY.etc.specparams.plotgroups,'together'), 1, 0);
    statval     = fastif(strcmpi(STUDY.etc.specparams.statistics,'param'), 1, 2);
    statcond    = fastif(strcmpi(STUDY.etc.specparams.statcond, 'on'), 1, 0);
    statgroup   = fastif(strcmpi(STUDY.etc.specparams.statgroup,'on'), 1, 0);
    
    uilist = { ...
        {'style' 'text'       'string' 'Frequency range (Hz)'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.specparams.freqrange) 'tag' 'freqrange' } ...
        {'style' 'text'       'string' 'Plot limit (uV)'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.specparams.ylim) 'tag' 'ylim' } ...
        {} {'style' 'checkbox'   'string' '' 'value' plotconditions 'enable' enablecond  'tag' 'plotconditions' } ...
        {'style' 'text'       'string' 'Plot conditions on the same panel' 'enable' enablecond } ...
        {} {'style' 'checkbox'   'string' '' 'value' plotgroups 'enable' enablegroup 'tag' 'plotgroups' } ...
        {'style' 'text'       'string' 'Plot groups on the same panel' 'enable' enablegroup } ...
        {} ...
        {'style' 'text'       'string' 'Statistics'} ...
        {'style' 'popupmenu'  'string' 'Parametric|Permutations' 'tag' 'statistics' 'value' statval 'listboxtop' statval } ...
        {'style' 'text'       'string' 'Threshold'} ...
        {'style' 'edit'       'string' threshstr 'tag' 'threshold' } ...
        {} {'style' 'checkbox'   'string' '' 'value' statcond  'enable' enablecond  'tag' 'statcond' } ...
        {'style' 'text'       'string' 'Compute condition statistics' 'enable' enablecond} ...
        {} {'style' 'checkbox'   'string' '' 'value' statgroup 'enable' enablegroup 'tag' 'statgroup' } ...
        {'style' 'text'       'string' 'Compute group statistics' 'enable' enablegroup } };
    
    geometry = { [ 1 1 1 1] [0.1 0.1 1] [0.1 0.1 1] [1] [1 1 1 1] [0.1 0.1 1] [0.1 0.1 1] };
    
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, ...
                                   'helpcom', 'pophelp(''std_specparams'')', ...
                                   'title', 'Set parameters for plotting specs -- pop_specparams()');

    if isempty(res), return; end;
    
    % decode inputs
    % -------------
    if res.plotgroups & res.plotconditions, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end;
    if res.statgroup, res.statgroup = 'on'; else res.statgroup = 'off'; end;
    if res.statcond , res.statcond  = 'on'; else res.statcond  = 'off'; end;
    if res.plotgroups, res.plotgroups = 'together'; else res.plotgroups = 'apart'; end;
    if res.plotconditions , res.plotconditions  = 'together'; else res.plotconditions  = 'apart'; end;
    res.freqrange = str2num( res.freqrange );
    res.ylim      = str2num( res.ylim );
    res.threshold = str2num( res.threshold );
    if isempty(res.threshold),res.threshold = NaN; end;
    if res.statistics == 1, res.statistics  = 'param'; 
    else                    res.statistics  = 'perm'; 
    end;
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( res.plotgroups, STUDY.etc.specparams.plotgroups), options = { options{:} 'plotgroups' res.plotgroups }; end;
    if ~strcmpi( res.plotconditions , STUDY.etc.specparams.plotconditions ), options = { options{:} 'plotconditions'  res.plotconditions  }; end;
    if ~strcmpi( res.statgroup, STUDY.etc.specparams.statgroup), options = { options{:} 'statgroup' res.statgroup }; end;
    if ~strcmpi( res.statcond , STUDY.etc.specparams.statcond ), options = { options{:} 'statcond'  res.statcond  }; end;
    if ~strcmpi( res.statistics, STUDY.etc.specparams.statistics ), options = { options{:} 'statistics' res.statistics }; end;
    if ~isequal(res.ylim, STUDY.etc.specparams.ylim),           options = { options{:} 'ylim' res.ylim      }; end;
    if ~isequal(res.freqrange, STUDY.etc.specparams.freqrange), options = { options{:} 'freqrange' res.freqrange }; end;
    if isnan(res.threshold) & ~isnan(STUDY.etc.specparams.threshold) | ...
            ~isnan(res.threshold) & isnan(STUDY.etc.specparams.threshold) | ...
                ~isnan(res.threshold) & res.threshold ~= STUDY.etc.specparams.threshold
                options = { options{:} 'threshold' res.threshold }; 
    end;
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
if ~isequal(STUDY.etc.specparams.freqrange, TMPSTUDY.etc.specparams.freqrange)
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
    if ~isfield(STUDY.etc.specparams, 'freqrange'),  STUDY.etc.specparams.freqrange = []; end;
    if ~isfield(STUDY.etc.specparams, 'ylim'     ),  STUDY.etc.specparams.ylim      = []; end;
    if ~isfield(STUDY.etc.specparams, 'statistics'), STUDY.etc.specparams.statistics = 'param'; end;
    if ~isfield(STUDY.etc.specparams, 'statgroup'),  STUDY.etc.specparams.statgroup = 'off'; end;
    if ~isfield(STUDY.etc.specparams, 'statcond' ),  STUDY.etc.specparams.statcond  = 'off'; end;
    if ~isfield(STUDY.etc.specparams, 'threshold' ), STUDY.etc.specparams.threshold = NaN; end;
    if ~isfield(STUDY.etc.specparams, 'plotgroups') , STUDY.etc.specparams.plotgroups = 'apart'; end;
    if ~isfield(STUDY.etc.specparams, 'plotconditions') ,  STUDY.etc.specparams.plotconditions  = 'apart'; end;
    if ~isfield(STUDY.etc.specparams, 'naccu')    ,  STUDY.etc.specparams.naccu     = []; end;
