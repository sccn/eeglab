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
%   'averagechan' - ['on'|'off'] average data channels when several are
%                  selected.
%
% See also: std_erpplot()
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

function [ STUDY, com ] = pop_erpparams(STUDY, varargin);

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    enablecond  = fastif(length(STUDY.design(STUDY.currentdesign).variable(1).value)>1, 'on', 'off');
    enablegroup = fastif(length(STUDY.design(STUDY.currentdesign).variable(2).value)>1, 'on', 'off');
    plotconditions    = fastif(strcmpi(STUDY.etc.erpparams.plotconditions, 'together'), 1, 0);
    plotgroups  = fastif(strcmpi(STUDY.etc.erpparams.plotgroups,'together'), 1, 0);
    radio_averagechan  = fastif(strcmpi(STUDY.etc.erpparams.averagechan,'on'), 1, 0);
    radio_scalptopo    = fastif(isempty(STUDY.etc.erpparams.topotime), 0, 1);
    if radio_scalptopo, radio_averagechan = 0; end;
    if radio_scalptopo+radio_averagechan == 0, radio_scalparray = 1; else radio_scalparray = 0; end;
        
    cb_radio = [ 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0);' ...
                 'set(gcbo, ''value'', 1);' ...
                 'set(findobj(gcbf, ''tag'', ''topotime''), ''string'', '''');' ];
    cb_edit  = [ 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0);' ...
                 'set(findobj(gcbf, ''tag'', ''scalptopotext''), ''value'', 1);' ];    
    uilist = { ...
        {'style' 'text'       'string' 'ERP plotting options' 'fontweight' 'bold' 'fontsize', 12} ...
        {} {'style' 'text'       'string' 'Frequency [low_Hz high_Hz]' } ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.timerange) 'tag' 'timerange' } ...
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
        {} {'style' 'radio'   'string' 'Plot channels in scalp array'    'value' radio_scalparray 'tag' 'scalparray'       'userdata' 'radio' 'callback' cb_radio} { } ...
        {} {'style' 'radio'   'string'  'Plot topography at time (ms)' 'value' radio_scalptopo  'tag' 'scalptopotext' 'userdata' 'radio' 'callback' cb_radio} ...
           {'style' 'edit'    'string' num2str(STUDY.etc.erpparams.topotime) 'tag' 'topotime' 'callback' cb_edit } ...
        {} {'style' 'radio'   'string' 'Average selected channels' 'value' radio_averagechan 'tag' 'averagechan' 'userdata' 'radio' 'callback' cb_radio} { } };
    cbline = [0.07 1.1];
    otherline = [ 0.07 0.6 .3];
    chanline  = [ 0.07 0.8 0.3];
    geometry = { 1 otherline otherline otherline 1 1 cbline cbline 1 1 chanline chanline chanline };
    geomvert = [1.2 1 1 1 0.5 1.2 1 1 0.5 1.2 1 1 1 ];
    
    % component plotting
    % ------------------
    if isnan(STUDY.etc.erpparams.topotime)
        geometry(end-4:end) = []; 
        geomvert(end-4:end) = []; 
        uilist(end-10:end) = [];
    end;
    
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'geomvert', geomvert, ...
                                            'title', 'ERP plotting options -- pop_erpparams()');
    if isempty(res), return; end;
    
    % decode inputs
    % -------------
    %if res.plotgroups & res.plotconditions, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end;
    if res.plotgroups, res.plotgroups = 'together'; else res.plotgroups = 'apart'; end;
    if res.plotconditions , res.plotconditions  = 'together'; else res.plotconditions  = 'apart'; end;
    if ~isfield(res, 'topotime'), res.topotime = STUDY.etc.erpparams.topotime;
    else res.topotime  = str2num( res.topotime );
    end;
    res.timerange = str2num( res.timerange );
    res.ylim      = str2num( res.ylim );
    res.filter    = str2num( res.filter );
    if ~isfield(res, 'averagechan'), res.averagechan = STUDY.etc.erpparams.averagechan;
    elseif res.averagechan, res.averagechan = 'on'; else res.averagechan = 'off';
    end;
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( char(res.filter), char(STUDY.etc.erpparams.filter)), options = { options{:} 'filter' res.filter }; end;
    if ~strcmpi( res.plotgroups, STUDY.etc.erpparams.plotgroups), options = { options{:} 'plotgroups' res.plotgroups }; end;
    if ~strcmpi( res.plotconditions , STUDY.etc.erpparams.plotconditions ), options = { options{:} 'plotconditions'  res.plotconditions  }; end;
    if ~isequal(res.ylim       , STUDY.etc.erpparams.ylim),      options = { options{:} 'ylim' res.ylim       }; end;
    if ~isequal(res.timerange  , STUDY.etc.erpparams.timerange)  , options = { options{:} 'timerange' res.timerange }; end;
    if ~isequal(res.averagechan, STUDY.etc.erpparams.averagechan), options = { options{:} 'averagechan' res.averagechan }; end;
    if (all(isnan(res.topotime)) & all(~isnan(STUDY.etc.erpparams.topotime))) | ...
            (all(~isnan(res.topotime)) & all(isnan(STUDY.etc.erpparams.topotime))) | ...
                (all(~isnan(res.topotime)) & ~isequal(res.topotime, STUDY.etc.erpparams.topotime))
        options = { options{:} 'topotime' res.topotime }; 
    end;
    if ~isempty(options)
        STUDY = pop_erpparams(STUDY, options{:});
        com = sprintf('STUDY = pop_erpparams(STUDY, %s);', vararg2str( options ));
    end;
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.erpparams), 'exact'))
                STUDY.etc.erpparams = setfield(STUDY.etc.erpparams, varargin{index}, varargin{index+1});
            end;
        end;
    end;
end;

% scan clusters and channels to remove erpdata info if timerange has changed
% ----------------------------------------------------------
if ~isequal(STUDY.etc.erpparams.timerange, TMPSTUDY.etc.erpparams.timerange)
    if isfield(STUDY.cluster, 'erpdata')
        for index = 1:length(STUDY.cluster)
            STUDY.cluster(index).erpdata  = [];
            STUDY.cluster(index).erptimes = [];
        end;
    end;
    if isfield(STUDY.changrp, 'erpdata')
        for index = 1:length(STUDY.changrp)
            STUDY.changrp(index).erpdata  = [];
            STUDY.changrp(index).erptimes = [];
        end;
    end;
end;

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'erpparams'), STUDY.etc.erpparams = []; end;
    if ~isfield(STUDY.etc.erpparams, 'topotime'),         STUDY.etc.erpparams.topotime = []; end;
    if ~isfield(STUDY.etc.erpparams, 'filter'),           STUDY.etc.erpparams.filter = []; end;
    if ~isfield(STUDY.etc.erpparams, 'timerange'),        STUDY.etc.erpparams.timerange = []; end;
    if ~isfield(STUDY.etc.erpparams, 'ylim'     ),        STUDY.etc.erpparams.ylim      = []; end;
    if ~isfield(STUDY.etc.erpparams, 'plotgroups') ,      STUDY.etc.erpparams.plotgroups = 'apart'; end;
    if ~isfield(STUDY.etc.erpparams, 'plotconditions') ,  STUDY.etc.erpparams.plotconditions  = 'apart'; end;
    if ~isfield(STUDY.etc.erpparams, 'averagechan') ,     STUDY.etc.erpparams.averagechan  = 'off'; end;

