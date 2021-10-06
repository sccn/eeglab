% pop_limo() - prepare and convert EEGLAB data and structure to be
%              processed by LIMO.
%
% Usage:
%   >> STUDY = pop_limo( STUDY, ALLEEG, 'key', val)
%   >> STUDY = pop_limo( STUDY, ALLEEG, 'dat', 'key', val)
%   >> STUDY = pop_limo( STUDY, ALLEEG, 'ica', 'key', val)
%
% Inputs:
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%   'dat'|'ica'  - show the interface for data channels or for ICA. The
%                  default is to use data.
%
% Optional inputs: same as std_limo()
%
% Graphical interface:
%   "See GLM factors" - [push button] See all the GLM factors or columns in
%                   the design matrix. This list depends on your study
%                   design. Adding contrast will change this list. In
%                   practice, it is important to understand that every
%                   single factor will be fitted by the GLM, and the choice
%                   of these factors is critical. It is recommended to
%                   include as many factors as possible, especially because
%                   this allows you to define new contrast later on without
%                   having to recompute the GLM on each subject (as long as
%                   these new contrast do not involve new factors).
%
%  "Interaction model for categorical indep. var." - When using more than 
%                   one categorical variable, clicking this option forces
%                   to have factors which are the conjunction of the
%                   different independent var. values. This is useful only 
%                   if you want to calculate interactions at the subject.
%                   The 'best' option is typically to have a design with
%                   all conditions (no factors) and create factors at the
%                   group level.
%
%  "Split regressions (continuous indep. var.)" - This options allows to
%                   split continuous variables for the different
%                   categorical variables. This is useful to compute 
%                   interaction between continuous and categorical
%                   variables. 
%
%  "Input data to use for the GLM" - [pop up menu] measure to use as input
%                   for the GLM. Currently, only "ERP" and "spectrum" are 
%                   supported.
%
%  "Optimization method" - [pop up menu] may be Ordinary Least Squares (OLS), 
%                   Weighted Least Squares (WLS), or Iterative Reweighted Least 
%                   Squares. WTS should be used as it is more robust to trials
%                   with different time course. OLS is a standard solution and 
%                   WLS/IRLS are solution that automatically weight down outliers
%                   (trials/data point).
%
%  "Options"        - [edit box] additional options. These are given
%                   directly as input to the std_limo.m function. They may
%                   be ['freqlim', value] frequency trimming in Hz or 
%                   ['timelim', value] time trimming in millisecond. These
%                   allow to speed up computation by not including all the
%                   data as input. Example: 'timelim', [-50 650] only
%                   include data within the -50 ms and 650 ms for ERPs.
%
% "Erase previous model" - [checkbox] erase previous model for this
%                   measure.
%
% Outputs:
%   STUDY       - an EEGLAB STUDY set of loaded EEG structures
%
% Author: Arnaud Delorme, SCCN, UCSD, 2015-
%         Cyril Pernet, LIMO Team - edit info and defaults
%
% See also: std_limo()

% Copyright (C) Arnaud Delorme
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

function [STUDY,com,limofiles] = pop_limo(STUDY, ALLEEG, measureflag, varargin)

if nargin < 2
    help pop_limo;
    return;
end

com = '';
if nargin > 2
    if strcmpi(measureflag,'dat')
        measureflagindx = 1;
    elseif strcmpi(measureflag,'ica')
        measureflagindx = 2;
    else
        varargin = { measureflag varargin{:} };
        measureflagindx = 1;
    end
else
    measureflagindx = 1;
end
    
if nargin < 4
    dataMeasures = { 'ERP' 'Spectrum' 'ERSP'};
    fileMeasures = { 'daterp' 'datspec' 'dattimef'; 'icaerp' 'icaspec' 'icatimef'};
    methods      = { 'WLS' 'OLS' 'IRLS'};
    cb_measure   = [ 'if get(gcbo, ''value'') == 1,' ...
                     '   set(findobj(gcbf, ''tag'', ''options''), ''string'', '''');' ...
                     'else,' ...
                     '   set(findobj(gcbf, ''tag'', ''options''), ''string'', ''''''freqlim'''', [1 25]'');' ...
                     'end;' ];
    cb_listfactors = [ 'pop_listfactors(STUDY, ''gui'', ''on'', ' ...
                               '''level'', ''one'',' ...
                               '''splitreg''   , fastif(get(findobj(gcbf, ''tag'', ''splitreg''   ), ''value''), ''on'', ''off''),' ...
                               '''interaction'', fastif(get(findobj(gcbf, ''tag'', ''interaction''), ''value''), ''on'', ''off''));' ];
    uilist = { ...
        {'style' 'text'       'string' 'LInear MOdeling of EEG data' 'fontweight' 'bold' 'fontsize', 12} ...
        {'style' 'pushbutton' 'string' 'See GLM variables' 'callback' cb_listfactors } ...
        {'style' 'checkbox'   'string' 'Interaction model for categorical indep. var.' 'value' 0 'tag' 'interaction' } ...
        {'style' 'checkbox'   'string' 'Split regressions (continuous indep. var.)' 'tag' 'splitreg' } {} ...
        {'style' 'text'       'string' 'Input data to use for the GLM' } ...
        {'style' 'popupmenu'  'string' dataMeasures 'tag' 'measure' 'callback' cb_measure} ...
        {'style' 'text'       'string' 'Optimization method' } ...
        {'style' 'popupmenu'  'string' methods 'tag' 'method' } ...
        {'style' 'text'       'string' 'Options' } ...
        {'style' 'edit'       'string' '' 'tag' 'options' } ...
        {'style' 'checkbox'   'string' 'Erase previous model' 'tag' 'erase' 'value' true } ...
        };
    
    cline = [1.1 0.8];
    geometry = { [1.6 1] 1 1 1 cline cline cline 1 };
    geomvert = [ 1 1 1     1 1 1     1     1 ];
        
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'geomvert', geomvert, ...
                                            'title', 'LInear MOdeling of EEG data -- pop_limo()', 'helpcom', 'pophelp(''pop_limo'');');
    if isempty(res), return; end
    opttmp  = eval( [ '{ ' res.options ' }' ]);
    if length(opttmp) > 0 && isnumeric(opttmp{1})
        error([ 'Wrong options. Options must be in the format' 10 '''key'', val. For example ''timelim'', [-100 600]' ]);
    end
    options = { 'method' methods{res.method} 'measure' fileMeasures{measureflagindx,res.measure} opttmp{:} ...
                'erase'       fastif(res.erase, 'on', 'off') ...
                'splitreg'    fastif(res.splitreg, 'on', 'off') ...
                'interaction' fastif(res.interaction, 'on', 'off') };
else
    options = varargin; 
end

[STUDY,limofiles] = std_limo(STUDY, ALLEEG, options{:});
com = sprintf('pop_limo(STUDY, ALLEEG, %s);', vararg2str(options));
