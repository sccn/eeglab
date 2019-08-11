% pop_erspparams() - Set plotting and statistics parameters for 
%                    computing and plotting STUDY mean (and optionally 
%                    single-trial) ERSP and ITC measures and measure 
%                    statistics. Settings are stored within the STUDY 
%                    structure (STUDY.etc.erspparams) which is used
%                    whenever plotting is performed by the function
%                    std_erspplot().
% Usage:    
%   >> STUDY = pop_erspparams(STUDY, 'key', 'val', ...);   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
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
%  'topofreq'    - [float] plot scalp map at specific frequencies. As above
%                  a frequency range may also be provided.
%  'subbaseline' - ['on'|'off'] subtract the same baseline across conditions 
%                  for ERSP (not ITC). When datasets with different conditions
%                  are recorded simultaneously, a common baseline spectrum 
%                  should be used. Note that this also affects the 
%                  results of statistics {default: 'on'}
%  'maskdata'    - ['on'|'off'] when threshold is not NaN, and 'groupstats'
%                  or 'condstats' (above) are 'off', masks the data 
%                  for significance. Deprecated.
%  'averagemode' - ['rms'|'ave'] average data channels when several are
%                  selected ('ave') or compute root mean square ('rms')
%                  which is the default.
%  'averagechan' - ['on'|'off'] average/rms data channels when several are
%                  selected ('on') or plot them individually ('off')
%
% See also: std_erspplot(), std_itcplot()
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2006-

% Copyright (C) Arnaud Delorme, 2006
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

function [ STUDY, com ] = pop_erspparams(STUDY, varargin);

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    subbaseline = fastif(strcmpi(STUDY.etc.erspparams.subbaseline,'on'), 1, 0);
    rmsFlag = fastif(strcmpi(STUDY.etc.erspparams.averagemode, 'rms'), 0, 1);
    icaFlag = fastif(isnan(STUDY.etc.erspparams.topotime), 1, 0);
    
    tmpTimeRange = STUDY.etc.erspparams.timerange;
    tmpFreqRange = STUDY.etc.erspparams.freqrange;
    if strcmpi(STUDY.etc.erspparams.averagechan,'off')
        if isempty(STUDY.etc.erspparams.topotime) || any(isnan(STUDY.etc.erspparams.topotime))
            multipleChansVal   = 1; % scalp array
        else
            multipleChansVal   = 2; % scalp topo
            tmpTimeRange = STUDY.etc.erspparams.topotime;
            tmpFreqRange = STUDY.etc.erspparams.topofreq;
        end
    else
        multipleChansVal   = 3; % average channels
    end
    
    cb_multiplechan    = [ 'if get(gcbo, ''value'') == 2 ' ...
                         '    if isempty(get(findobj(gcbf, ''tag'', ''timerange''), ''string'')) || isempty(get(findobj(gcbf, ''tag'', ''freqrange''), ''string'')),' ...
                         '       set(gcbo, ''value'', 1);' ...
                         '       warndlg2([ ''Select time and frequency range first for plotting'' 10 ''topographies, then select that setting again.'' ]);' ...
                         '    end;' ...
                         'else,' ...
                         '    if ~isempty(get(findobj(gcbf, ''tag'', ''timerange''), ''string'')) && length(unique(str2num(get(findobj(gcbf, ''tag'', ''timerange''), ''string'')))) ==1,' ...
                         '       set(findobj(gcbf, ''tag'', ''timerange''), ''string'', '''');' ...
                         '    end;' ...
                         '    if ~isempty(get(findobj(gcbf, ''tag'', ''freqrange''), ''string'')) && length(unique(str2num(get(findobj(gcbf, ''tag'', ''freqrange''), ''string'')))) ==1,' ...
                         '       set(findobj(gcbf, ''tag'', ''freqrange''), ''string'', '''');' ...
                         '    end;' ...
                         'end;' ];    
    uilist = { ...
        {'style' 'text'       'string' 'ERSP/ITC plotting options' 'fontweight' 'bold' 'tag', 'ersp' } ...
        {'style' 'text'       'string' 'Time range in ms [Low High]'} ...
        {'style' 'edit'       'string' num2str(tmpTimeRange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Freq. range in Hz [Low High]'} ...
        {'style' 'edit'       'string' num2str(tmpFreqRange) 'tag' 'freqrange' } ...
        {'style' 'text'       'string' 'Power limits in dB [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erspparams.ersplim) 'tag' 'ersplim' } ...
        {'style' 'text'       'string' 'ITC limit (0-1) [High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erspparams.itclim) 'tag' 'itclim' } ...
        {} {'style' 'checkbox'   'string' 'Average time/freq images instead of RMS' 'value' rmsFlag 'tag' 'averagemode' }  ...
        {} {'style' 'checkbox'   'string' 'Common ERSP baseline across factors' 'value' subbaseline 'tag' 'subbaseline' }  ...
        {} ...
        {'style' 'text'       'string' 'Multiple channels selection' 'fontweight' 'bold' 'tag', 'spec' 'fontsize', 12} ...
        {} {'style' 'popupmenu'  'string' { 'Plot channels individually' 'Plot averaged topography over time/freq range' 'Average/RMS of selected channels' } 'value' multipleChansVal 'tag' 'multiplechan' 'callback' cb_multiplechan } {} };
    evalstr = 'set(findobj(gcf, ''tag'', ''ersp''), ''fontsize'', 12);';
    cbline = [0.07 1.1];
    otherline = [ 0.6 .4 ];
    chanline  = [ 0.07 0.8];
    geometry = { 1 otherline otherline otherline otherline cbline cbline 1 1 chanline 1 };
    
    if icaFlag
        uilist = uilist(1:end-5);
        geometry = geometry(1:end-4);
    end
    
    [~, ~, ~, res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'skipline', 'off', ...
                                            'title', 'Set ERSP/ITC plotting parameters -- pop_erspparams()', 'eval', evalstr);
    if isempty(res), return; end
    
    % decode input
    % ------------
    if res.subbaseline, res.subbaseline = 'on'; else res.subbaseline = 'off'; end
    if res.averagemode, res.averagemode = 'ave'; else res.averagemode = 'rms'; end
    if ~isfield(res, 'multiplechan') res.multiplechan = 0; end
    res.timerange = str2num( res.timerange );
    res.freqrange = str2num( res.freqrange );
    res.ersplim   = str2num( res.ersplim );
    res.itclim    = str2num( res.itclim );
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( res.averagemode , STUDY.etc.erspparams.averagemode ), options = { options{:} 'averagemode' res.averagemode }; end
    if ~strcmpi( res.subbaseline , STUDY.etc.erspparams.subbaseline ), options = { options{:} 'subbaseline' res.subbaseline }; end
    if ~isequal(res.ersplim  , STUDY.etc.erspparams.ersplim),   options = { options{:} 'ersplim'   res.ersplim   }; end
    if ~isequal(res.itclim   , STUDY.etc.erspparams.itclim),    options = { options{:} 'itclim'    res.itclim    }; end
    if ~isequal(res.timerange, STUDY.etc.erspparams.timerange) && res.multiplechan ~= 2, options = { options{:} 'timerange' res.timerange }; end
    if ~isequal(res.freqrange, STUDY.etc.erspparams.freqrange) && res.multiplechan ~= 2, options = { options{:} 'freqrange' res.freqrange }; end
    
    % mutliple channel option
    % -----------------------
    if res.multiplechan == 1
        if ~isequal('off', STUDY.etc.erspparams.averagechan), options = { options{:} 'averagechan' 'off' }; end
        if ~isempty(       STUDY.etc.erspparams.topotime),    options = { options{:} 'topotime' [] }; end
        if ~isempty(       STUDY.etc.erspparams.topofreq),    options = { options{:} 'topofreq' [] }; end
    elseif res.multiplechan == 2
        if ~isequal('off', STUDY.etc.erspparams.averagechan), options = { options{:} 'averagechan' 'off' }; end
        if ~isequal(res.timerange, STUDY.etc.erspparams.topotime) options = { options{:} 'topotime' res.timerange }; end
        if ~isequal(res.freqrange, STUDY.etc.erspparams.topofreq) options = { options{:} 'topofreq' res.freqrange }; end
        if ~isequal([], STUDY.etc.erspparams.timerange) options = { options{:} 'timerange' [] }; end
        if ~isequal([], STUDY.etc.erspparams.freqrange) options = { options{:} 'freqrange' [] }; end
    elseif res.multiplechan > 2
        if ~isempty(       STUDY.etc.erspparams.topotime),    options = { options{:} 'topotime' [] }; end
        if ~isempty(       STUDY.etc.erspparams.topofreq),    options = { options{:} 'topofreq' [] }; end
        if ~isequal('on', STUDY.etc.erspparams.averagechan), options = { options{:} 'averagechan' 'on' }; end
    end
    
    % execute options
    % ---------------
    if ~isempty(options)
        STUDY = pop_erspparams(STUDY, options{:});
        com = sprintf('STUDY = pop_erspparams(STUDY, %s);', vararg2str( options ));
    end
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.erspparams), 'exact'))
                STUDY.etc.erspparams = setfield(STUDY.etc.erspparams, varargin{index}, varargin{index+1});
            end
        end
    end
end

% scan clusters and channels to remove erspdata info if timerange etc. have changed
% ---------------------------------------------------------------------------------
if ~isequal(STUDY.etc.erspparams.timerange, TMPSTUDY.etc.erspparams.timerange) || ... 
    ~isequal(STUDY.etc.erspparams.freqrange, TMPSTUDY.etc.erspparams.freqrange) || ... 
    ~isequal(STUDY.etc.erspparams.subbaseline, TMPSTUDY.etc.erspparams.subbaseline)
    rmfields = { 'erspdata' 'ersptimes' 'erspfreqs' 'erspbase' 'erspdatatrials' 'ersptimes' 'erspfreqs' 'erspsubjinds' 'ersptrialinfo' ...
                 'itcdata'  'itctimes'  'itcfreqs'             'itcdatatrials'  'itctimes'  'itcfreqs'  'itcsubjinds'  'itctrialinfo' };
    for iField = 1:length(rmfields)
        if isfield(STUDY.cluster, rmfields{iField})
            STUDY.cluster = rmfield(STUDY.cluster, rmfields{iField});
        end
        if isfield(STUDY.changrp, rmfields{iField})
            STUDY.changrp = rmfield(STUDY.changrp, rmfields{iField});
        end
    end
end

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'erspparams'), STUDY.etc.erspparams = []; end
    if ~isfield(STUDY.etc.erspparams, 'topotime'),     STUDY.etc.erspparams.topotime = []; end
    if ~isfield(STUDY.etc.erspparams, 'topofreq'),     STUDY.etc.erspparams.topofreq = []; end
    if ~isfield(STUDY.etc.erspparams, 'timerange'),    STUDY.etc.erspparams.timerange = []; end
    if ~isfield(STUDY.etc.erspparams, 'freqrange'),    STUDY.etc.erspparams.freqrange = []; end
    if ~isfield(STUDY.etc.erspparams, 'ersplim' ),     STUDY.etc.erspparams.ersplim   = []; end
    if ~isfield(STUDY.etc.erspparams, 'itclim' ),      STUDY.etc.erspparams.itclim    = []; end
    if ~isfield(STUDY.etc.erspparams, 'maskdata' ),    STUDY.etc.erspparams.maskdata  = 'off'; end %deprecated
    if ~isfield(STUDY.etc.erspparams, 'averagemode' ), STUDY.etc.erspparams.averagemode  = 'rms'; end
    if ~isfield(STUDY.etc.erspparams, 'averagechan' ), STUDY.etc.erspparams.averagechan  = 'off'; end
    if ~isfield(STUDY.etc.erspparams, 'subbaseline' ),  STUDY.etc.erspparams.subbaseline = 'off'; end

