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
%                   groups {default: 'off'}
%  'condstats'    - ['on'|'off'] Compute statistics across data
%                   conditions {default: 'off'}
%  'statistics'   - ['param'|'perm'|'bootstrap'] Type of statistics to compute
%                   'param' for parametric (t-test/anova); 'perm' for
%                   permutation-based and 'bootstrap' for bootstrap
%                   {default: 'param'}
%  'singletrials' - ['on'|'off'] use single trials to compute statistics.
%                   This requires the measure to be computed with the
%                   'savetrials', 'on' option.
%  'mode'         - ['eeglab'|'fieldtrip'] use either EEGLAB or Fieldtrip
%                   statistics.
%  'effect'       - ['main'|'marginal'] compute main effect or marginal
%                   effects.
%
% EEGLAB statistics:
%  'method'      - ['param'|'perm'|'bootstrap'] statistical 
%                  method. See help statcond for more information.
%  'naccu'       - [integer] Number of surrogate data averages to use in
%                  surrogate statistics. For instance, if p<0.01,
%                  use naccu>200. For p<0.001, naccu>2000. If a 'threshold'
%                  (not NaN) is set below and 'naccu' is too low, it will
%                  be automatically increased. (This keyword is currently
%                  only modifiable from the command line, not from the gui).
%  'alpha'       - [NaN|alpha] Significance threshold (0<alpha<<1). Value
%                  NaN will plot p-values for each time and/or frequency
%                  on a different axis. If alpha is used, significant time
%                  and/or frequency regions will be indicated either on
%                  a separate axis or (whenever possible) along with the
%                  data {default: NaN}
%  'mcorrect'    - ['fdr'|'holms'|'bonferoni'|'none'] correction for multiple
%                  comparisons. 'fdr' uses false discovery rate. See the fdr 
%                  function for more information. Default is
%                  none'.
%
% Fieldtrip statistics:
%  'fieldtripmethod' - ['analytic'|'montecarlo'] statistical 
%                  method. See help statcond for more information.
%  'fieldtripnaccu' - [integer] Number of surrogate data averages to use in
%                  surrogate statistics.
%  'fieldtripalpha' - [alpha] Significance threshold (0<alpha<<1). This
%                  parameter is mandatory. Default is 0.05.
%  'fieldtripmcorrect' - ['cluster'|'max'|'fdr'|'holms'|'bonferoni'|'none'] 
%                  correction for multiple comparisons. See help 
%                  ft_statistics_montecarlo for more information.
%  'fieldtripclusterparam - [string] parameters for clustering. See help 
%                  ft_statistics_montecarlo for more information. 
%  'fieldtripchannelneighbor - [struct] channel neighbor structure.
%  'fieldtripchannelneighborparam' - [string] parameters for channel 
%                  neighbor. See help ft_statistics_montecarlo for more 
%                  information.
%
% Legacy parameters:
%  'threshold'   - now 'alpha'
%  'statistics'  - now 'method'
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2010-

% Copyright (C) Arnaud Delorme, 2010
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

function [STUDY, com] = pop_statparams(STUDY, varargin);

com = '';
if isfield(STUDY, 'etc')
    if ~isfield(STUDY.etc, 'statistics') STUDY.etc.statistics = default_stats([]);
    else                                 STUDY.etc.statistics = default_stats(STUDY.etc.statistics);
    end
    if length(varargin) == 1 && strcmpi(varargin{1}, 'default')
        return;
    end
end

if isempty(varargin) && ~isempty(STUDY)
    if ~exist('ft_freqstatistics'), fieldtripInstalled = false; else fieldtripInstalled = true; end
    opt.enablesingletrials = 'on';
    
    % encode parameters
    % -----------------
    paramstruct         = STUDY.etc.statistics;
    eeglabStatvalues    = { 'param' 'perm' 'bootstrap' };
    fieldtripStatvalues = { 'analytic' 'montecarlo' };
    mCorrectList        = { 'none' 'bonferoni' 'holms' 'fdr' 'max' 'cluster' };
    marginal   = fastif(strcmpi(paramstruct.effect, 'main'), 0, 1);
    condstats   = fastif(strcmpi(paramstruct.condstats, 'on'), 1, 0);
    groupstats  = fastif(strcmpi(paramstruct.groupstats,'on'), 1, 0);
    statmode    = fastif(strcmpi(paramstruct.singletrials,'on'), 1, 0);
    eeglabThresh    = fastif(isnan(paramstruct.eeglab.alpha),'exact', num2str(paramstruct.eeglab.alpha));
    fieldtripThresh = fastif(isnan(paramstruct.fieldtrip.alpha),'exact', num2str(paramstruct.fieldtrip.alpha));
    eeglabStat      = strmatch(paramstruct.eeglab.method, eeglabStatvalues);
    fieldtripStat   = strmatch(paramstruct.fieldtrip.method, fieldtripStatvalues);
    if isempty(eeglabStat)   , error('Unknown statistical method for EEGLAB'); end
    if isempty(fieldtripStat), error('Unknown statistical method for Fieldtrip'); end
    eeglabRand      = fastif(isempty(paramstruct.eeglab.naccu), 'auto', num2str(paramstruct.eeglab.naccu));
    fieldtripRand   = fastif(isempty(paramstruct.fieldtrip.naccu), 'auto', num2str(paramstruct.fieldtrip.naccu));
    eeglabMcorrect    = strmatch(paramstruct.eeglab.mcorrect, mCorrectList);
    fieldtripMcorrect = strmatch(paramstruct.fieldtrip.mcorrect, mCorrectList);
    fieldtripClust    = paramstruct.fieldtrip.clusterparam;
    fieldtripChan     = paramstruct.fieldtrip.channelneighborparam;
    combootstrap =  [ 'warndlg2( strvcat(''Bootstrap is selected. Bootstrap is provided for legacy reasons.'',' ...
        '''It overestimates significance by a factor of 2 for paired data '',' ...
        '''(unpaired data is working properly). It is advised to use'',' ...
        '''permutation instead of bootstrap if you have paired data.'') );' ];
    cb_bootstrap     = [ 'if get(gcbo, ''value'') == 3,' combootstrap ' end;' ];
    cb_help_neighbor = 'pophelp(''ft_prepare_neighbours'');';
    cb_help_cluster  = 'pophelp(''ft_statistics_montecarlo'');';
    cb_textSyntax    = 'try, eval( [ ''{'' get(gcbo, ''string'') ''};'' ]); catch, warndlg2(''Syntax error''); end;';
    
    % callback for randomization selection
    % ------------------------------------
    cbSelectRandEeglab   = [ 'set(findobj(gcbf, ''tag'', ''eeglabnaccu''),    ''enable'', ''on'');' ...
        'set(findobj(gcbf, ''tag'', ''eeglabnaccutext''),''enable'', ''on'');' ];
    cbUnselectRandEeglab = [ 'set(findobj(gcbf, ''tag'', ''eeglabnaccu''),    ''enable'', ''off'');' ...
        'set(findobj(gcbf, ''tag'', ''eeglabnaccutext''),''enable'', ''off'');' ];
    cbSelectRandFieldtrip  = [ 'set(findobj(gcbf, ''tag'', ''fieldtripnaccu''),    ''enable'', ''on'');' ...
        'set(findobj(gcbf, ''tag'', ''fieldtripnaccutext''),''enable'', ''on'');' ];
    cbUnselectRandFieldtrip = [ 'set(findobj(gcbf, ''tag'', ''fieldtripnaccu''),    ''enable'', ''off'');' ...
        'set(findobj(gcbf, ''tag'', ''fieldtripnaccutext''),''enable'', ''off'');' ];
    cbSetFullMcorrectFieldtrip = 'set(findobj(gcbf, ''tag'', ''fieldtripmcorrect''), ''string'', ''Do not correct for multiple comparisons|Use Bonferoni correction|Use Holms correction|Use FDR correction|Use max correction|Use cluster correction (CC)'');';
    cbUnsetFullMcorrectFieldtrip = 'set(findobj(gcbf, ''tag'', ''fieldtripmcorrect''), ''value'', min(4, get(findobj(gcbf, ''tag'', ''fieldtripmcorrect''), ''value'')), ''string'', ''Do not correct for multiple comparisons|Use Bonferoni correction|Use Holms correction|Use FDR correction'');';
    cb_eeglab_statlist    = [ 'if get(findobj(gcbf, ''tag'', ''stateeglab'' ),''value'') > 1,' cbSelectRandEeglab ',else,' cbUnselectRandEeglab ',end;' ];
    cb_fieldtrip_statlist = [ 'if get(findobj(gcbf, ''tag'', ''statfieldtrip'' ),''value'') > 1,' cbSelectRandFieldtrip cbSetFullMcorrectFieldtrip ',else,' cbUnselectRandFieldtrip cbUnsetFullMcorrectFieldtrip ',end;' ];
    
    % callback for activating clusters inputs
    % ---------------------------------------
    cb_select_cluster   = [ 'set(findobj(gcbf, ''tag'', ''clustertext1''),''enable'', ''on'');' ...
        'set(findobj(gcbf, ''tag'', ''clustertext2''),''enable'', ''on'');' ...
        'set(findobj(gcbf, ''tag'', ''clusterhelp1''),''enable'', ''on'');' ...
        'set(findobj(gcbf, ''tag'', ''clusterhelp2''),''enable'', ''on'');' ...
        'set(findobj(gcbf, ''tag'', ''clusterchan'' ),''enable'', ''on'');' ...
        'set(findobj(gcbf, ''tag'', ''clusterstat'' ),''enable'', ''on'');' ];
    cb_unselect_cluster = [ 'set(findobj(gcbf, ''tag'', ''clustertext1''),''enable'', ''off'');' ...
        'set(findobj(gcbf, ''tag'', ''clustertext2''),''enable'', ''off'');' ...
        'set(findobj(gcbf, ''tag'', ''clusterhelp1''),''enable'', ''off'');' ...
        'set(findobj(gcbf, ''tag'', ''clusterhelp2''),''enable'', ''off'');' ...
        'set(findobj(gcbf, ''tag'', ''clusterchan'' ),''enable'', ''off'');' ...
        'set(findobj(gcbf, ''tag'', ''clusterstat'' ),''enable'', ''off'');' ];
    cb_fieldtrip_mcorrect = [ cb_fieldtrip_statlist 'if get(findobj(gcbf, ''tag'', ''fieldtripmcorrect'' ),''value'') == 6,' cb_select_cluster ',else,' cb_unselect_cluster ',end;' cb_fieldtrip_statlist]; % cb_fieldtrip_statlist repeated on purpose
    
    % callback for activating eeglab/fieldtrip
    % ----------------------------------------
    enable_eeglab = [  'set(findobj(gcbf, ''userdata'', ''eeglab'')   ,''enable'', ''on'');' ...
                       'set(findobj(gcbf, ''userdata'', ''fieldtrip''),''enable'', ''off'');' ...
                       'set(findobj(gcbf, ''tag'', ''but_eeglab'')   ,''value'', 1);' ...
                       'set(findobj(gcbf, ''tag'', ''but_fieldtrip''),''value'', 0);' cb_eeglab_statlist ];
    
    enable_fieldtrip=[ 'if get(findobj(gcbf, ''tag'', ''condstats''), ''value'') && get(findobj(gcbf, ''tag'', ''groupstats''), ''value''),' ...
        enable_eeglab ...
        'warndlg2(strvcat(''Switching to EEGLAB statistics since'',''Fieldtrip cannot perform 2-way statistics''));' ...
        'else,' ...
        '  set(findobj(gcbf, ''userdata'', ''eeglab'')   ,''enable'', ''off'');' ...
        '  set(findobj(gcbf, ''userdata'', ''fieldtrip''),''enable'', ''on'');' ...
        '  set(findobj(gcbf, ''tag'', ''but_fieldtrip''),''value'', 1);' ...
        '  set(findobj(gcbf, ''tag'', ''but_eeglab'')   ,''value'', 0);' cb_fieldtrip_mcorrect ...
        'end;'];
    
    cb_select_fieldtrip = [ 'if get(findobj(gcbf, ''tag'', ''but_fieldtrip''),''value''),' enable_fieldtrip ',else,' enable_eeglab    ',end;' ];
    cb_select_eeglab    = [ 'if get(findobj(gcbf, ''tag'', ''but_eeglab''),''value''),,'   enable_eeglab    ',else,' enable_fieldtrip ',end;' ];
    
    if strcmpi(paramstruct.mode, 'eeglab'), evalstr = enable_eeglab;
    else                                    evalstr = enable_fieldtrip;
    end
    inds = findstr('gcbf', evalstr);
    evalstr(inds+2) = [];
    
    % special case if Fieldtrip is not installed
    if ~fieldtripInstalled
        strFieldtrip = 'Use Fieldtrip statistics (to use install "Fieldtrip-lite" using File > Manage EEGLAB extensions)';
        fieldtripEnable = 'off';
        cb_select_eeglab = 'set(findobj(gcbf, ''tag'', ''but_eeglab''),''value'', 1)';
    else
        strFieldtrip = 'Use Fieldtrip statistics';
        fieldtripEnable = 'on';
    end
    
    opt.uilist = { ...
        {'style' 'text'        'string' 'General statistical parameters' 'fontweight' 'bold' } ...
        {} {'style' 'checkbox' 'string' 'Compute 1st independent variable statistics if any' 'value' condstats  'callback' cb_select_fieldtrip 'tag' 'condstats' } ...
        {} {'style' 'checkbox' 'string' 'Compute 2nd independent variable statistics if any' 'value' groupstats 'callback' cb_select_fieldtrip 'tag' 'groupstats' } ...
        {} {'style' 'checkbox' 'string' 'Plot marginal statistics (unckeck to plot main effect)' 'value' marginal  'tag' 'marginal' } ...
        {} {'style' 'checkbox' 'string' 'Use single trials for statistics (not recommended if more than 1 subject in the study)' 'value' statmode 'tag' 'singletrials' 'enable' opt.enablesingletrials } ...
        {} ...
        {'style' 'checkbox'     'string' 'Use EEGLAB statistics' 'fontweight' 'bold' 'tag' 'but_eeglab' 'callback' cb_select_eeglab } ...
        {} {'style' 'popupmenu' 'string' { 'Use parametric statistics' 'Use permutation statistics' 'Use bootstrap statistics' } 'tag' 'stateeglab' 'value' eeglabStat 'listboxtop' eeglabStat 'callback' [ cb_eeglab_statlist cb_bootstrap ] 'userdata' 'eeglab' } ...
        {'style' 'text'         'string' 'Statistical threshold (p-value)' 'userdata' 'eeglab'} ...
        {'style' 'edit'         'string' eeglabThresh 'tag' 'eeglabalpha' 'userdata' 'eeglab' } ...
        {} {'style' 'popupmenu' 'string' { 'Do not correct for multiple comparisons' 'Use Bonferoni correction' 'Use Holms correction' 'Use FDR correction' } 'value' eeglabMcorrect 'listboxtop' eeglabMcorrect 'tag' 'eeglabmcorrect' 'userdata' 'eeglab'  } ...
        {'style' 'text' 'string' '   Randomization (n)' 'userdata' 'eeglab' 'tag' 'eeglabnaccutext' } ...
        {'style' 'edit' 'string' eeglabRand          'userdata' 'eeglab' 'tag' 'eeglabnaccu' } ...
        {} ...
        {'style' 'checkbox'     'string' strFieldtrip 'enable' fieldtripEnable 'fontweight' 'bold' 'tag' 'but_fieldtrip' 'callback' cb_select_fieldtrip } ...
        {} {'style' 'popupmenu' 'string' { 'Use analytic/parametric statistics' 'Use montecarlo/permutation statistics' } 'tag' 'statfieldtrip' 'value' fieldtripStat 'listboxtop' fieldtripStat 'callback' cb_fieldtrip_mcorrect 'userdata' 'fieldtrip' } ...
        {'style' 'text'        'string' 'Statistical threshold (p-value)' 'userdata' 'fieldtrip' } ...
        {'style' 'edit'        'string' fieldtripThresh 'tag' 'fieldtripalpha' 'userdata' 'fieldtrip' } ...
        {} {'style' 'popupmenu' 'string' { 'Do not correct for multiple comparisons' 'Use Bonferoni correction' 'Use Holms correction' 'Use FDR correction' 'Use max correction' 'Use cluster correction (CC)' } 'tag' 'fieldtripmcorrect' 'value' fieldtripMcorrect 'listboxtop' fieldtripMcorrect 'callback' cb_fieldtrip_mcorrect 'userdata' 'fieldtrip' } ...
        {'style' 'text' 'string' '   Randomization (n)' 'userdata' 'fieldtrip' 'tag' 'fieldtripnaccutext' } ...
        {'style' 'edit' 'string' fieldtripRand       'userdata' 'fieldtrip' 'tag' 'fieldtripnaccu' } ...
        {} {'style' 'text' 'string' 'CC channel neighbor parameters'       'userdata' 'fieldtrip' 'tag' 'clustertext1' } ...
        { 'style' 'edit' 'string' fieldtripChan           'userdata' 'fieldtrip' 'tag' 'clusterchan' 'callback' cb_textSyntax } ...
        { 'style' 'pushbutton' 'string' 'help' 'callback' cb_help_neighbor 'userdata' 'fieldtrip' 'tag' 'clusterhelp1' } ...
        {} {'style' 'text' 'string' 'CC clustering parameters'             'userdata' 'fieldtrip' 'tag' 'clustertext2' }  ...
        { 'style' 'edit' 'string' fieldtripClust        'userdata' 'fieldtrip' 'tag' 'clusterstat' 'callback' cb_textSyntax } ...
        { 'style' 'pushbutton' 'string' 'help' 'callback' cb_help_cluster  'userdata' 'fieldtrip' 'tag' 'clusterhelp2' } ...
        };
    
    if eeglabStat == 3,
        eval(combootstrap);
    end
    
    cbline = [0.07 1.1];
    otherline = [ 0.7 0.6 .5];
    eeglabline = [ 0.7 0.6 .5];
    opt.geometry = { [1] cbline cbline cbline cbline [1] [1] [0.07 0.51 0.34 0.13] [0.07 0.6 0.25 0.13] ...
        [1] [1] [0.07 0.51 0.34 0.13] [0.07 0.6 0.25 0.13] [0.13 0.4 0.4 0.1] [0.13 0.4 0.4 0.1]  };
    
    [out_param, userdat, tmp, res] = inputgui( 'geometry' , opt.geometry, 'uilist', opt.uilist, ...
                                            'title', 'Set statistical parameters -- pop_statparams()','eval', evalstr);
    if isempty(res), return; end
    
    % decode parameters
    % ----------------
    if res.marginal,    res.effect       = 'marginal'; else res.effect = 'main';  end
    if res.groupstats,   res.groupstats   = 'on';  else res.groupstats = 'off';  end
    if res.condstats ,   res.condstats    = 'on';  else res.condstats  = 'off';  end
    if res.singletrials, res.singletrials = 'on';  else res.singletrials = 'off';  end
    res.eeglabalpha    = str2num(res.eeglabalpha);
    res.fieldtripalpha = str2num(res.fieldtripalpha);
    if isempty(res.eeglabalpha)   ,res.eeglabalpha = NaN;    end
    if isempty(res.fieldtripalpha),res.fieldtripalpha = NaN; end
    res.stateeglab         = eeglabStatvalues{res.stateeglab};
    res.statfieldtrip      = fieldtripStatvalues{res.statfieldtrip};
    res.eeglabmcorrect     = mCorrectList{res.eeglabmcorrect};
    res.fieldtripmcorrect  = mCorrectList{res.fieldtripmcorrect};
    res.mode               = fastif(res.but_eeglab, 'eeglab', 'fieldtrip');
    res.eeglabnaccu        = str2num(res.eeglabnaccu);
    if ~ischar(res.fieldtripnaccu) || ~strcmpi(res.fieldtripnaccu, 'all')
        res.fieldtripnaccu     = str2num(res.fieldtripnaccu);
    end
    
    % build command call
    % ------------------
    options = {};
    if strcmp(res.stateeglab, 'param' ) && exist('fcdf') ~= 2
        fprintf(['statcond(): EEGLAB parametric testing requires fcdf() \n' ...
            '                 from the Matlab Statstical Toolbox. Running\n' ...
            '                 nonparametric permutation tests instead.\n']);
        res.stateeglab = 'perm';
    end
    if ~strcmpi( res.effect,       paramstruct.effect),        options = { options{:} 'effect' res.effect }; end
    if ~strcmpi( res.groupstats,   paramstruct.groupstats),    options = { options{:} 'groupstats' res.groupstats }; end
    if ~strcmpi( res.condstats ,   paramstruct.condstats ),    options = { options{:} 'condstats'  res.condstats  }; end
    if ~strcmpi( res.singletrials, paramstruct.singletrials ), options = { options{:} 'singletrials'  res.singletrials }; end
    if ~strcmpi( res.mode              , paramstruct.mode),                options = { options{:} 'mode' res.mode }; end % statistics
    if ~isequal( res.eeglabnaccu       , paramstruct.eeglab.naccu),        options = { options{:} 'naccu'    res.eeglabnaccu }; end
    if ~strcmpi( res.stateeglab        , paramstruct.eeglab.method),       options = { options{:} 'method' res.stateeglab }; end % statistics
    if ~strcmpi( res.eeglabmcorrect    , paramstruct.eeglab.mcorrect),     options = { options{:} 'mcorrect' res.eeglabmcorrect }; end
    if ~isequal( res.fieldtripnaccu    , paramstruct.fieldtrip.naccu),     options = { options{:} 'fieldtripnaccu' res.fieldtripnaccu }; end
    if ~strcmpi( res.statfieldtrip     , paramstruct.fieldtrip.method),    options = { options{:} 'fieldtripmethod' res.statfieldtrip }; end
    if ~strcmpi( res.fieldtripmcorrect , paramstruct.fieldtrip.mcorrect),  options = { options{:} 'fieldtripmcorrect' res.fieldtripmcorrect }; end
    if ~strcmpi( res.clusterstat       , paramstruct.fieldtrip.clusterparam),  options = { options{:} 'fieldtripclusterparam' res.clusterstat }; end
    if ~strcmpi( res.clusterchan       , paramstruct.fieldtrip.channelneighborparam),  options = { options{:} 'fieldtripchannelneighborparam' res.clusterchan }; end
    if ~(isnan(res.eeglabalpha(1)) && isnan(paramstruct.eeglab.alpha(1))) && ~isequal(res.eeglabalpha, paramstruct.eeglab.alpha) % threshold
        options = { options{:} 'alpha' res.eeglabalpha };
    end
    if ~(isnan(res.fieldtripalpha(1)) && isnan(paramstruct.fieldtrip.alpha(1))) && ~isequal(res.fieldtripalpha, paramstruct.fieldtrip.alpha) % threshold
        options = { options{:} 'fieldtripalpha' res.fieldtripalpha };
    end
    
    if ~isempty(options)
        STUDY = pop_statparams(STUDY, options{:});
        com = sprintf('STUDY = pop_statparams(STUDY, %s);', vararg2str( options ));
    end
else
    % interpret parameters
    % --------------------
    if isfield(STUDY, 'etc')
         paramstruct = STUDY.etc.statistics; isstudy = true;
    else paramstruct = STUDY;
         if isempty(paramstruct), paramstruct = default_stats([]); end
         isstudy = false;
    end
    
    if isempty(varargin) || strcmpi(varargin{1}, 'default')
        paramstruct = default_stats(paramstruct);
    else
        for index = 1:2:length(varargin)
            v = varargin{index};
            if strcmpi(v, 'statistics'), v = 'method'; end % backward compatibility
            if strcmpi(v, 'threshold' ), v = 'alpha';  end % backward compatibility
            
            if strcmpi(v, 'alpha') || strcmpi(v, 'method') || strcmpi(v, 'naccu') || strcmpi(v, 'mcorrect') 
                paramstruct = setfield(paramstruct, 'eeglab', v, varargin{index+1});
            elseif ~isempty(findstr('fieldtrip', v))
                v2 = v(10:end);
                paramstruct = setfield(paramstruct, 'fieldtrip', v2, varargin{index+1});
                if strcmpi(v2, 'channelneighborparam')
                    paramstruct.fieldtrip.channelneighbor = []; % reset neighbor matrix if parameter change
                end
            else
                if (~isempty(paramstruct) && ~isempty(strmatch(v, fieldnames(paramstruct), 'exact'))) || ~isstudy
                    paramstruct = setfield(paramstruct, v, varargin{index+1});
                end
            end
        end
    end
    
    if isfield(STUDY, 'etc')
         STUDY.etc.statistics = paramstruct; 
    else STUDY = paramstruct;
    end
end

% default parameters
% ------------------
function paramstruct = default_stats(paramstruct)

if ~isfield(paramstruct, 'effect'),        paramstruct.effect     = 'main'; end
if ~isfield(paramstruct, 'groupstats'),    paramstruct.groupstats = 'off';  end
if ~isfield(paramstruct, 'condstats' ),    paramstruct.condstats  = 'off';  end
if ~isfield(paramstruct, 'singletrials' ), paramstruct.singletrials = 'off'; end
if ~isfield(paramstruct, 'mode' ),         paramstruct.mode         = 'eeglab'; end
if ~isfield(paramstruct, 'eeglab'),        paramstruct.eeglab       = []; end
if ~isfield(paramstruct, 'fieldtrip'),     paramstruct.fieldtrip    = []; end
if ~isfield(paramstruct.eeglab, 'naccu'),    paramstruct.eeglab.naccu = []; end
if ~isfield(paramstruct.eeglab, 'alpha' ),   paramstruct.eeglab.alpha = NaN; end
if ~isfield(paramstruct.eeglab, 'method'),   paramstruct.eeglab.method = 'param'; end
if ~isfield(paramstruct.eeglab, 'mcorrect'), paramstruct.eeglab.mcorrect = 'none'; end
if ~isfield(paramstruct.fieldtrip, 'naccu'),  paramstruct.fieldtrip.naccu = []; end
if ~isfield(paramstruct.fieldtrip, 'method'), paramstruct.fieldtrip.method = 'analytic'; end
if ~isfield(paramstruct.fieldtrip, 'alpha'),  paramstruct.fieldtrip.alpha  = NaN; end
if ~isfield(paramstruct.fieldtrip, 'mcorrect'),   paramstruct.fieldtrip.mcorrect = 'none'; end
if ~isfield(paramstruct.fieldtrip, 'clusterparam'),   paramstruct.fieldtrip.clusterparam = '''clusterstatistic'',''maxsum'''; end
if ~isfield(paramstruct.fieldtrip, 'channelneighbor'),   paramstruct.fieldtrip.channelneighbor = []; end
if ~isfield(paramstruct.fieldtrip, 'channelneighborparam'),   paramstruct.fieldtrip.channelneighborparam = '''method'',''triangulation'''; end
if strcmpi(paramstruct.eeglab.mcorrect, 'benferoni'), paramstruct.eeglab.mcorrect = 'bonferoni'; end
if strcmpi(paramstruct.eeglab.mcorrect,    'no'), paramstruct.eeglab.mcorrect = 'none'; end
if strcmpi(paramstruct.fieldtrip.mcorrect, 'no'), paramstruct.fieldtrip.mcorrect = 'none'; end
