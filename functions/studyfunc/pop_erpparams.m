% pop_erpparams() - Set plotting and statistics parameters for cluster ERP 
%                   plotting
% Usage:    
%   >> STUDY = pop_erpparams(STUDY, 'key', 'val');   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% Input:
%   'topotime'   - [real] Plot ERP scalp maps at one specific latency (ms).
%                   A latency range [min max] may also be defined (the 
%                   ERP is then averaged over the interval) {default: []}
%   'filter'     - [real] low pass filter the ERP curves at a given 
%                  frequency threshold. Default is no filtering.
%   'timerange'  - [min max] ERP plotting latency range in ms. 
%                  {default: the whole epoch}
%   'ylim'       - [min max] ERP limits in microvolts {default: from data}
%   'plotgroups' - ['together'|'apart'] 'together' -> plot subject groups 
%                  on the same axis in different colors, else ('apart') 
%                  on different axes. {default: 'apart'}
%   'plotconditions' - ['together'|'apart'] 'together' -> plot conditions 
%                  on the same axis in different colors, else ('apart') 
%                  on different axes. Note: Keywords 'plotgroups' and 
%                  'plotconditions' may not both be set to 'together'. 
%                  {default: 'together'}
%   'averagechan' - ['rms'|'on'|'off'] average data channels when several are
%                  selected ('on') or compute root mean square ('rms').
%
% See also: std_erpplot()
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

function [ STUDY, com ] = pop_erpparams(STUDY, varargin)

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    enablecond  = 'off';
    enablegroup = 'off';
    if length(STUDY.design(STUDY.currentdesign).variable) > 0 && length(STUDY.design(STUDY.currentdesign).variable(1).value)>1, enablecond  = 'on'; end
    if length(STUDY.design(STUDY.currentdesign).variable) > 1 && length(STUDY.design(STUDY.currentdesign).variable(2).value)>1, enablegroup = 'on'; end;   
    plotconditions     = fastif(strcmpi(STUDY.etc.erpparams.plotconditions, 'together'), 1, 0);
    plotgroups         = fastif(strcmpi(STUDY.etc.erpparams.plotgroups,'together'), 1, 0);
    tmpTimeRange = STUDY.etc.erpparams.timerange;
    if strcmpi(STUDY.etc.erpparams.averagechan,'off')
        if isempty(STUDY.etc.erpparams.topotime) || any(isnan(STUDY.etc.erpparams.topotime))
            multipleChansVal   = 1; % scalp array
        else
            multipleChansVal   = 2; % scalp topo
            tmpTimeRange = STUDY.etc.erpparams.topotime;
        end
    else
        if strcmpi(STUDY.etc.erpparams.averagechan,'on')
            multipleChansVal   = 3; % average channels
        else
            multipleChansVal   = 4; % root mean square
        end
    end
        
    cb_radio = [ 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0);' ...
                 'set(gcbo, ''value'', 1);' ...
                 'set(findobj(gcbf, ''tag'', ''topotime''), ''string'', '''');' ];
    cb_edit  = [ 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0);' ...
                 'set(findobj(gcbf, ''tag'', ''scalptopotext''), ''value'', 1);' ];    
    cb_multiplechan    = [ 'if get(gcbo, ''value'') == 2 ' ...
                         '    if isempty(get(findobj(gcbf, ''tag'', ''timerange''), ''string'')),' ...
                         '       set(gcbo, ''value'', 1);' ...
                         '       warndlg2([ ''Select time range first for plotting topography,'' 10 ''then select that setting again.'' ]);' ...
                         '    end;' ...
                         'else,' ...
                         '    if ~isempty(get(findobj(gcbf, ''tag'', ''timerange''), ''string'')) && length(unique(str2num(get(findobj(gcbf, ''tag'', ''timerange''), ''string'')))) ==1,' ...
                         '       set(findobj(gcbf, ''tag'', ''timerange''), ''string'', '''');' ...
                         '    end;' ...
                         'end;' ];
    uilist = { ...
           {'style' 'text'       'string' 'ERP plotting options' 'fontweight' 'bold' 'fontsize', 12} ...
        {} {'style' 'text'       'string' 'Time range (ms) [low high]' } ...
           {'style' 'edit'       'string' num2str(tmpTimeRange) 'tag' 'timerange' } ...
        {} {'style' 'text'       'string' 'Plot limits [low high]'} ...
           {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.ylim) 'tag' 'ylim' } ...
        {} {'style' 'text'       'string' 'Lowpass plotted data [Hz]' } ...
           {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.filter) 'tag' 'filter' } ...
        {} ...
           {'style' 'text'       'string' 'ERP plotting format' 'fontweight' 'bold' 'fontsize', 12} ...
        {} {'style' 'checkbox'   'string' 'Plot first variable on the same panel' 'value' plotconditions 'enable' enablecond  'tag' 'plotconditions' } ...
        {} {'style' 'checkbox'   'string' 'Plot second variable on the same panel' 'value' plotgroups 'enable' enablegroup 'tag' 'plotgroups' } ...
        {} ...
           {'style' 'text'       'string' 'Multiple channels selection' 'fontweight' 'bold' 'tag', 'spec' 'fontsize', 12} ...
        {} {'style' 'popupmenu'  'string' { 'Plot channels individually' 'Plot averaged topography over time range' 'Average potential of selected channels' 'Compute RMS of selected channels' } 'value' multipleChansVal 'tag' 'multiplechan' 'callback' cb_multiplechan } };
    cbline = [0.07 1.1];
    otherline = [ 0.07 0.6 .3];
    chanline  = [ 0.07 0.8];
    geometry = { 1 otherline otherline otherline 1 1 cbline cbline 1 1 chanline };
    geomvert = [1.2 1 1 1 0.5 1.2 1 1 0.5 1.2 1 ];
    
    % component plotting
    % ------------------
    if isnan(STUDY.etc.erpparams.topotime)
        geometry(end-2:end) = []; 
        geomvert(end-2:end) = []; 
        uilist(end-3:end) = [];
    end
    
    [out_param, userdat, tmp, res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'geomvert', geomvert, ...
                                            'title', 'ERP plotting options -- pop_erpparams()');
    if isempty(res), return; end
    if ~isfield(res, 'multiplechan'), res.multiplechan = 0; end
    
    % decode inputs
    % -------------
    %if res.plotgroups && res.plotconditions, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end
    if res.plotgroups, res.plotgroups = 'together'; else res.plotgroups = 'apart'; end
    if res.plotconditions , res.plotconditions  = 'together'; else res.plotconditions  = 'apart'; end
    res.timerange = str2num( res.timerange );
    res.ylim      = str2num( res.ylim );
    res.filter    = str2num( res.filter );
    if ~isempty(diff(res.timerange))
        if diff(res.timerange) < 5
            fprintf(2, 'Time range is less than 5 ms; are you sure you entered the time range in milliseconds, not seconds?\n');
        end
    end

    % build command call
    % ------------------
    options = {};
    if ~strcmpi( char(res.filter), char(STUDY.etc.erpparams.filter)), options = { options{:} 'filter' res.filter }; end
    if ~strcmpi( res.plotgroups, STUDY.etc.erpparams.plotgroups), options = { options{:} 'plotgroups' res.plotgroups }; end
    if ~strcmpi( res.plotconditions , STUDY.etc.erpparams.plotconditions ), options = { options{:} 'plotconditions'  res.plotconditions  }; end
    if ~isequal(res.ylim       , STUDY.etc.erpparams.ylim),      options = { options{:} 'ylim' res.ylim       }; end
    if ~isequal(res.timerange  , STUDY.etc.erpparams.timerange) &&  res.multiplechan ~= 2, options = { options{:} 'timerange' res.timerange }; end
    
    % mutliple channel option
    % -----------------------
    if res.multiplechan == 1
        if ~isequal('off', STUDY.etc.erpparams.averagechan), options = { options{:} 'averagechan' 'off' }; end
        if ~isempty(       STUDY.etc.erpparams.topotime),    options = { options{:} 'topotime' [] }; end
    elseif res.multiplechan == 2
        if ~isequal('off', STUDY.etc.erpparams.averagechan), options = { options{:} 'averagechan' 'off' }; end
        if ~isequal(res.timerange, STUDY.etc.erpparams.topotime) options = { options{:} 'topotime' res.timerange }; end
        if ~isequal([], STUDY.etc.erpparams.timerange) options = { options{:} 'timerange' [] }; end
        if isempty(res.timerange)
            disp('Warning: you must select a time range to plot scalp topographies, plotting individual channels instead');
        end
    elseif res.multiplechan > 2
        if ~isempty(       STUDY.etc.erpparams.topotime),    options = { options{:} 'topotime' [] }; end
        if res.multiplechan == 3
            if ~isequal('on', STUDY.etc.erpparams.averagechan), options = { options{:} 'averagechan' 'on' }; end
        else
            if ~isequal('rms', STUDY.etc.erpparams.averagechan), options = { options{:} 'averagechan' 'rms' }; end
        end
    end
    
    % execute option
    % --------------
    if ~isempty(options)
        STUDY = pop_erpparams(STUDY, options{:});
        com = sprintf('STUDY = pop_erpparams(STUDY, %s);', vararg2str( options ));
    end
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.erpparams), 'exact'))
                STUDY.etc.erpparams = setfield(STUDY.etc.erpparams, varargin{index}, varargin{index+1});
            end
        end
    end
end

% scan clusters and channels to remove erpdata info if timerange has changed
% ----------------------------------------------------------
if ~isequal(STUDY.etc.erpparams.timerange, TMPSTUDY.etc.erpparams.timerange)
    rmfields = { 'erpdata' 'erptimes' };
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
    if ~isfield(STUDY.etc, 'erpparams'), STUDY.etc.erpparams = []; end
    if ~isfield(STUDY.etc.erpparams, 'topotime'),         STUDY.etc.erpparams.topotime = []; end
    if ~isfield(STUDY.etc.erpparams, 'filter'),           STUDY.etc.erpparams.filter = []; end
    if ~isfield(STUDY.etc.erpparams, 'timerange'),        STUDY.etc.erpparams.timerange = []; end
    if ~isfield(STUDY.etc.erpparams, 'ylim'     ),        STUDY.etc.erpparams.ylim      = []; end
    if ~isfield(STUDY.etc.erpparams, 'plotgroups') ,      STUDY.etc.erpparams.plotgroups = 'apart'; end
    if ~isfield(STUDY.etc.erpparams, 'plotconditions') ,  STUDY.etc.erpparams.plotconditions  = 'apart'; end
    if ~isfield(STUDY.etc.erpparams, 'averagechan') ,     STUDY.etc.erpparams.averagechan  = 'off'; end
    if ~isfield(STUDY.etc.erpparams, 'detachplots') ,     STUDY.etc.erpparams.detachplots  = 'on'; end % deprecated

