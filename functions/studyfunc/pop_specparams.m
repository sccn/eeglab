% pop_specparams() - Set plotting and statistics parameters for computing
%                    STUDY component spectra.
% Usage:    
%       >> STUDY = pop_specparams(STUDY, 'key', 'val');   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% Plot options:
%   'topofreq'   - [real] Plot Spectrum scalp maps at one specific freq. (Hz).
%                  A frequency range [min max] may also be defined (the 
%                  spectrum is then averaged over the interval) {default: []}
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
%   'subtractsubjectmean' - ['on'|'off'] subtract individual subject mean
%                  from each spectrum before plotting and computing
%                  statistics. Default is 'off'.
%   'averagechan' - ['rms'|'on'|'off'] average data channels when several are
%                  selected ('on') or compute root mean square ('rms').
%
% See also: std_specplot()
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2006-

% Copyright (C) Arnaud Delorme, CERCO
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

function [ STUDY, com ] = pop_specparams(STUDY, varargin);

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    enablecond  = 'off';
    enablegroup = 'off';
    if length(STUDY.design(STUDY.currentdesign).variable) > 0 && length(STUDY.design(STUDY.currentdesign).variable(1).value)>1, enablecond  = 'on'; end
    if length(STUDY.design(STUDY.currentdesign).variable) > 1 && length(STUDY.design(STUDY.currentdesign).variable(2).value)>1, enablegroup = 'on'; end;   
    plotconditions     = fastif(strcmpi(STUDY.etc.specparams.plotconditions, 'together'), 1, 0);
    plotgroups         = fastif(strcmpi(STUDY.etc.specparams.plotgroups,'together'), 1, 0);
    submean            = fastif(strcmpi(STUDY.etc.specparams.subtractsubjectmean,'on'), 1, 0);
    tmpFreqRange = STUDY.etc.specparams.freqrange;
    if strcmpi(STUDY.etc.specparams.averagechan,'off')
        if isempty(STUDY.etc.specparams.topofreq) || any(isnan(STUDY.etc.specparams.topofreq))
            multipleChansVal   = 1; % scalp array
        else
            multipleChansVal   = 2; % scalp topo
            tmpFreqRange = STUDY.etc.specparams.topofreq;
        end
    else
        if strcmpi(STUDY.etc.specparams.averagechan,'on')
            multipleChansVal   = 3; % average channels
        else
            multipleChansVal   = 4; % root mean square
        end
    end
            
    cb_radio = [ 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0);' ...
                 'set(gcbo, ''value'', 1);' ...
                 'set(findobj(gcbf, ''tag'', ''topofreq''), ''string'', '''');' ];
    cb_edit  = [ 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0);' ...
                 'set(findobj(gcbf, ''tag'', ''scalptopotext''), ''value'', 1);' ];
    cb_multiplechan    = [ 'if get(gcbo, ''value'') == 2 ' ...
                         '    if isempty(get(findobj(gcbf, ''tag'', ''freqrange''), ''string'')),' ...
                         '       set(gcbo, ''value'', 1);' ...
                         '       warndlg2([ ''Select frequency range first for plotting topography,'' 10 ''then select that setting again.'' ]);' ...
                         '    end;' ...
                         'else,' ...
                         '    if ~isempty(get(findobj(gcbf, ''tag'', ''freqrange''), ''string'')) && length(unique(str2num(get(findobj(gcbf, ''tag'', ''freqrange''), ''string'')))) ==1,' ...
                         '       set(findobj(gcbf, ''tag'', ''freqrange''), ''string'', '''');' ...
                         '    end;' ...
                         'end;' ];
    
    uilist = { ...
           {'style' 'text'       'string' 'Spectrum plotting options' 'fontweight' 'bold' 'fontsize', 12} ...
        {} {'style' 'text'       'string' 'Frequency [low_Hz high_Hz]' } ...
           {'style' 'edit'       'string' num2str(tmpFreqRange) 'tag' 'freqrange' } ...
        {} {'style' 'text'       'string' 'Plot limits [low high]'} ...
           {'style' 'edit'       'string' num2str(STUDY.etc.specparams.ylim) 'tag' 'ylim' } ...
        {} {'style' 'checkbox'   'string' 'Subtract individual subject mean spectrum' 'value' submean 'tag' 'submean' } ...
        {} ...
           {'style' 'text'       'string' 'Spectrum plotting format' 'fontweight' 'bold' 'fontsize', 12} ...
        {} {'style' 'checkbox'   'string' 'Plot first variable on the same panel' 'value' plotconditions 'enable' enablecond  'tag' 'plotconditions' } ...
        {} {'style' 'checkbox'   'string' 'Plot second variable on the same panel' 'value' plotgroups 'enable' enablegroup 'tag' 'plotgroups' } ...
        {} ...
           {'style' 'text'       'string' 'Multiple channels selection' 'fontweight' 'bold' 'tag', 'spec' 'fontsize', 12} ...
        {} {'style' 'popupmenu'  'string' { 'Plot channels individually' 'Plot averaged topography over frequency range' 'Average power of selected channels' 'Compute RMS power of selected channels' } 'value' multipleChansVal 'tag' 'multiplechan' 'callback' cb_multiplechan } };
    cbline = [0.07 1.1];
    otherline = [ 0.07 0.6 .3];
    chanline  = [ 0.07 0.8];
    geometry = { 1 otherline otherline cbline 1 1 cbline cbline 1 1 chanline };
    geomvert = [1.2 1 1 1 0.5 1.2 1 1 0.5 1.2 1 ];
    
    % component plotting
    % ------------------
    if isnan(STUDY.etc.specparams.topofreq)
        geometry(end-2:end) = []; 
        geomvert(end-2:end) = []; 
        uilist(end-3:end) = [];
    end
    
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'geomvert', geomvert, ...
                                            'title', 'Spectrum plotting options -- pop_specparams()');
    if isempty(res), return; end
    if ~isfield(res, 'multiplechan'), res.multiplechan = 0; end
    
    % decode inputs
    % -------------
    %if res.plotgroups && res.plotconditions, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end
    if res.submean   , res.submean    = 'on'; else res.submean    = 'off'; end
    if res.plotgroups, res.plotgroups = 'together'; else res.plotgroups = 'apart'; end
    if res.plotconditions , res.plotconditions  = 'together'; else res.plotconditions  = 'apart'; end
    res.freqrange = str2num( res.freqrange );
    res.ylim      = str2num( res.ylim );
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( res.plotgroups, STUDY.etc.specparams.plotgroups), options = { options{:} 'plotgroups' res.plotgroups }; end
    if ~strcmpi( res.plotconditions , STUDY.etc.specparams.plotconditions ), options = { options{:} 'plotconditions'  res.plotconditions  }; end
    if ~strcmpi( res.submean   , STUDY.etc.specparams.subtractsubjectmean ), options = { options{:} 'subtractsubjectmean'  res.submean  }; end
    if ~isequal(res.ylim, STUDY.etc.specparams.ylim),               options = { options{:} 'ylim' res.ylim      }; end
    if ~isequal(res.freqrange, STUDY.etc.specparams.freqrange) && res.multiplechan ~= 2,     options = { options{:} 'freqrange' res.freqrange }; end
    
    % mutliple channel option
    % -----------------------
    if res.multiplechan == 1
        if ~isequal('off', STUDY.etc.specparams.averagechan), options = { options{:} 'averagechan' 'off' }; end
        if ~isempty(       STUDY.etc.specparams.topofreq),    options = { options{:} 'topofreq' [] }; end
    elseif res.multiplechan == 2
        if ~isequal('off', STUDY.etc.specparams.averagechan), options = { options{:} 'averagechan' 'off' }; end
        if ~isequal(res.freqrange, STUDY.etc.specparams.topofreq) options = { options{:} 'topofreq' res.freqrange }; end
        if ~isequal([], STUDY.etc.specparams.freqrange) options = { options{:} 'freqrange' [] }; end
        if isempty(res.freqrange)
            disp('Warning: you must select a frequency range to plot scalp topographies, plotting individual channels instead');
        end
    elseif res.multiplechan > 2
        if ~isempty(       STUDY.etc.specparams.topofreq),    options = { options{:} 'topofreq' [] }; end
        if res.multiplechan == 3
            if ~isequal('on', STUDY.etc.specparams.averagechan), options = { options{:} 'averagechan' 'on' }; end
        else
            if ~isequal('rms', STUDY.etc.specparams.averagechan), options = { options{:} 'averagechan' 'rms' }; end
        end
    end
    
    % execute option
    % --------------
    if ~isempty(options)
        STUDY = pop_specparams(STUDY, options{:});
        com = sprintf('STUDY = pop_specparams(STUDY, %s);', vararg2str( options ));
    end
else
    
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.specparams), 'exact'))
                STUDY.etc.specparams = setfield(STUDY.etc.specparams, varargin{index}, varargin{index+1});
            end
        end
    end
end

% scan clusters and channels to remove specdata info if freqrange has changed
% ----------------------------------------------------------
if ~isequal(STUDY.etc.specparams.freqrange, TMPSTUDY.etc.specparams.freqrange) || ...
    ~isequal(STUDY.etc.specparams.subtractsubjectmean, TMPSTUDY.etc.specparams.subtractsubjectmean)
    rmfields = { 'specdata' 'specfreqs' };
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
    if ~isfield(STUDY.etc, 'specparams'), STUDY.etc.specparams = []; end
    if ~isfield(STUDY.etc.specparams, 'topofreq'),             STUDY.etc.specparams.topofreq = []; end
    if ~isfield(STUDY.etc.specparams, 'freqrange'),            STUDY.etc.specparams.freqrange = []; end
    if ~isfield(STUDY.etc.specparams, 'ylim'     ),            STUDY.etc.specparams.ylim      = []; end
    if ~isfield(STUDY.etc.specparams, 'subtractsubjectmean' ), STUDY.etc.specparams.subtractsubjectmean  = 'off'; end
    if ~isfield(STUDY.etc.specparams, 'plotgroups'),           STUDY.etc.specparams.plotgroups = 'apart'; end
    if ~isfield(STUDY.etc.specparams, 'plotconditions'),       STUDY.etc.specparams.plotconditions  = 'apart'; end
    if ~isfield(STUDY.etc.specparams, 'averagechan') ,         STUDY.etc.specparams.averagechan  = 'off'; end
    if ~isfield(STUDY.etc.specparams, 'detachplots') ,         STUDY.etc.specparams.detachplots  = 'on'; end % deprecated
