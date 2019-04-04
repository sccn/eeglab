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
%                  for significance.
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
    vis = fastif(isnan(STUDY.etc.erspparams.topotime), 'off', 'on');
    
    uilist = { ...
        {'style' 'text'       'string' 'ERSP/ITC plotting options' 'fontweight' 'bold' 'tag', 'ersp' } ...
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
        {} {'style' 'checkbox'   'string' 'Compute common ERSP baseline for all conditions/factors' 'value' subbaseline 'tag' 'subbaseline' }  };
    evalstr = 'set(findobj(gcf, ''tag'', ''ersp''), ''fontsize'', 12);';
    cbline = [0.07 1.1];
    otherline = [ 0.6 .4 0.6 .4];
    geometry = { 1 otherline otherline otherline cbline };
    enablecond  = 'off';
    enablegroup = 'off';
    if length(STUDY.design(STUDY.currentdesign).variable) > 0 && length(STUDY.design(STUDY.currentdesign).variable(1).value)>1, enablecond  = 'on'; end
    if length(STUDY.design(STUDY.currentdesign).variable) > 1 && length(STUDY.design(STUDY.currentdesign).variable(2).value)>1, enablegroup = 'on'; end;   
    
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'skipline', 'off', ...
                                            'title', 'Set ERSP/ITC plotting parameters -- pop_erspparams()', 'eval', evalstr);
    if isempty(res), return; end
    
    % decode input
    % ------------
    if res.subbaseline, res.subbaseline = 'on'; else res.subbaseline = 'off'; end
    res.topotime  = str2num( res.topotime );
    res.topofreq  = str2num( res.topofreq );
    res.timerange = str2num( res.timerange );
    res.freqrange = str2num( res.freqrange );
    res.ersplim   = str2num( res.ersplim );
    res.itclim    = str2num( res.itclim );
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( res.subbaseline , STUDY.etc.erspparams.subbaseline ), options = { options{:} 'subbaseline' res.subbaseline }; end
    if ~isequal(res.topotime , STUDY.etc.erspparams.topotime),  options = { options{:} 'topotime'   res.topotime  }; end
    if ~isequal(res.topofreq , STUDY.etc.erspparams.topofreq),  options = { options{:} 'topofreq'   res.topofreq  }; end
    if ~isequal(res.ersplim  , STUDY.etc.erspparams.ersplim),   options = { options{:} 'ersplim'   res.ersplim   }; end
    if ~isequal(res.itclim   , STUDY.etc.erspparams.itclim),    options = { options{:} 'itclim'    res.itclim    }; end
    if ~isequal(res.timerange, STUDY.etc.erspparams.timerange), options = { options{:} 'timerange' res.timerange }; end
    if ~isequal(res.freqrange, STUDY.etc.erspparams.freqrange), options = { options{:} 'freqrange' res.freqrange }; end
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
    if ~isfield(STUDY.etc.erspparams, 'maskdata' ),    STUDY.etc.erspparams.maskdata  = 'off'; end; %deprecated
    if ~isfield(STUDY.etc.erspparams, 'subbaseline' ),  STUDY.etc.erspparams.subbaseline = 'off'; end

