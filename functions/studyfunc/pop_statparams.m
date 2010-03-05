% pop_statparams() - helper function for pop_erspparams, pop_erpparams, and
%                    pop_specparams.
%
% Usage:    
%   >> struct = pop_statparams(struct, 'default');   
%   >> struct = pop_statparams(struct, 'key', 'val', ...);   
%
% Inputs:
%   struct        - parameter structure. When called with the 'default'
%                   option only, the function creates all the fields in 
%                   the structure and populate these fields with default 
%                   values.
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
%  'singletrials' - ['on'|'off'] use single trials to compute statistics.
%                  This requires the measure to be computed with the
%                  'savetrials', 'on' option.
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2010-

% Copyright (C) Arnaud Delorme, 2010
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
% Revision 1.2  2010/02/24 15:18:38  claire
% fix singletrials
%
% Revision 1.1  2010/02/24 10:52:36  arno
% Implemented new single trial statistics
%

function [paramstruct, res, options] = pop_statparams(paramstruct, varargin);

paramstruct = default_stats(paramstruct);
options = {};
if length(varargin) == 0 || strcmpi(varargin{1}, 'default')
    return;
end;

[opt ignore] = finputcheck(varargin, { 'uilist'   'cell'    []    {};
                                       'enablegroup'   'string' { 'on' 'off' }  'on';
                                       'enablecond'    'string' { 'on' 'off' }  'off';
                                       'geometry' 'cell'    []    {} }, 'pop_params', 'ignore');
if isstr(opt), error(opt); end;

condstats   = fastif(strcmpi(paramstruct.condstats, 'on'), 1, 0);
mcorrect    = fastif(strcmpi(paramstruct.mcorrect,  'fdr'), 1, 0);
groupstats  = fastif(strcmpi(paramstruct.groupstats,'on'), 1, 0);
statmode    = fastif(strcmpi(paramstruct.singletrials,'on'), 1, 0);
threshstr   = fastif(isnan(paramstruct.threshold),'', num2str(paramstruct.threshold));
if strcmpi(paramstruct.statistics,'param'),    statval = 1;
elseif strcmpi(paramstruct.statistics,'perm'), statval = 2;
else                                           statval = 3;
end;
cb_maskdata = [ 'tmpcond  = get(findobj(gcbf, ''tag'', ''condstats'') , ''value'');' ...
                'tmpgroup = get(findobj(gcbf, ''tag'', ''groupstats''), ''value'');' ...
                'tmpplot  = get(findobj(gcbf, ''tag'', ''maskdata'') , ''value'');' ...
                'if tmpcond & tmpgroup & tmpplot,' ...
                '    warndlg2(strvcat(''Cannot mask time/freq. image if both statistics for conditions'',' ...
                '           ''and statistics for groups are used.''));' ...
                '    set(gcbo, ''value'', 0);' ...
                'end;' ...
                'clear tmpcond tmpgroup tmpplot;' ];
                
opt.uilist = { opt.uilist{:} ...
        {} ...
        {'style' 'text'        'string' 'Statistical method to use'} ...
        {'style' 'popupmenu'   'string' 'Parametric|Permutation|Bootstrap' 'tag' 'statistics' 'value' statval 'listboxtop' statval } ...
        {'style' 'text'        'string' 'Statistical threshold (p<)' } ...
        {'style' 'edit'        'string' threshstr 'tag' 'threshold' } ...
        {} {'style' 'checkbox' 'string' 'Compute condition statistics' 'value' condstats  'enable' opt.enablecond  'tag' 'condstats' } ...
        {} {'style' 'checkbox' 'string' 'Compute group statistics' 'value' groupstats 'enable' opt.enablegroup 'tag' 'groupstats' } ...
        {} {'style' 'checkbox' 'string' 'Use single trials (when available)' 'value' statmode 'tag' 'singletrials' } ...
        {} {'style' 'checkbox' 'string' 'Use False Discovery Rate to correct for multiple comparisons' 'value' mcorrect 'tag' 'mcorrect' } };

cbline = [0.07 1.1];
otherline = [ 0.7 .5 0.6 .5];
opt.geometry = { opt.geometry{:} [1] otherline cbline cbline cbline cbline };
 
[out_param userdat tmp res] = inputgui( 'geometry' , opt.geometry, 'uilist', opt.uilist, ignore{:});
if isempty(res), return; end;

if res.groupstats,   res.groupstats   = 'on';  else res.groupstats = 'off';  end;
if res.condstats ,   res.condstats    = 'on';  else res.condstats  = 'off';  end;
if res.mcorrect,     res.mcorrect     = 'fdr'; else res.mcorrect   = 'none'; end;
if res.singletrials, res.singletrials = 'on';  else res.singletrials   = 'off';  end;
res.threshold = str2num(res.threshold);
if isempty(res.threshold),res.threshold = NaN; end;
if res.statistics == 1, res.statistics  = 'param'; 
elseif res.statistics == 2, res.statistics  = 'perm'; 
else res.statistics  = 'bootstrap'; 
end;

% build command call
% ------------------
if ~strcmpi( res.groupstats,   paramstruct.groupstats),    options = { options{:} 'groupstats' res.groupstats }; end;
if ~strcmpi( res.condstats ,   paramstruct.condstats ),    options = { options{:} 'condstats'  res.condstats  }; end;
if ~strcmpi( res.singletrials, paramstruct.singletrials ), options = { options{:} 'singletrials'  res.singletrials }; end;
if ~strcmpi( res.statistics,   paramstruct.statistics ),   options = { options{:} 'statistics' res.statistics }; end;
if ~strcmpi( res.mcorrect,     paramstruct.mcorrect),      options = { options{:} 'mcorrect' res.mcorrect }; end;
if (~isempty(res.threshold) && isnan(res.threshold(1)) && ~isnan(paramstruct.threshold(1))) || ...
        (~isempty(res.threshold) && ~isnan(res.threshold(1)) && isnan(paramstruct.threshold(1))) || ...
            (~isempty(res.threshold) && ~isnan(res.threshold(1)) && ~isequal(res.threshold, paramstruct.threshold))
            options = { options{:} 'threshold' res.threshold }; 
end;

function paramstruct = default_stats(paramstruct)
    if ~isfield(paramstruct, 'statistics'),    paramstruct.statistics = 'param'; end;
    if ~isfield(paramstruct, 'groupstats'),    paramstruct.groupstats = 'off';  end;
    if ~isfield(paramstruct, 'condstats' ),    paramstruct.condstats  = 'off';  end;
    if ~isfield(paramstruct, 'subbaseline' ),  paramstruct.subbaseline = 'on';  end;
    if ~isfield(paramstruct, 'singletrials' ), paramstruct.singletrials = 'off'; end;
    if ~isfield(paramstruct, 'threshold' ),    paramstruct.threshold = NaN; end;
    if ~isfield(paramstruct, 'mcorrect' ),     paramstruct.mcorrect   = 'none'; end;    
    if ~isfield(paramstruct, 'naccu')    ,     paramstruct.naccu     = []; end;

