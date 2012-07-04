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
    vis = fastif(isnan(STUDY.etc.erpparams.topotime), 'off', 'on');
    
    uilist = { ...
        {'style' 'text'       'string' 'ERP plotting options' 'fontweight' 'bold' 'tag', 'erp' } ...
        {'style' 'text'       'string' 'Time range in ms [low high]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.timerange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Plot limits in uV [low high]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.ylim) 'tag' 'ylim' } ...
        {'style' 'text'       'string' 'Plot scalp map at latency [ms]' 'enable' vis } ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.topotime) 'tag' 'topotime' 'enable' vis } ...
        {'style' 'text'       'string' 'Display filter in Hz [high]' } ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.filter) 'tag' 'filter' } ...
        {} {'style' 'checkbox'   'string' 'Plot first variable on the same panel' 'value' plotconditions 'enable' enablecond  'tag' 'plotconditions' } ...
        {} {'style' 'checkbox'   'string' 'Plot second variable on the same panel' 'value' plotgroups 'enable' enablegroup 'tag' 'plotgroups' } };
    evalstr = 'set(findobj(gcf, ''tag'', ''erp''), ''fontsize'', 12);';
    cbline = [0.07 1.1];
    otherline = [ 0.7 .5 0.6 .5];
    geometry = { 1 otherline otherline cbline cbline };
    
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'skipline', 'off', ...
                                            'title', 'Set ERP plotting parameters -- pop_erpparams()', 'eval', evalstr);
    if isempty(res), return; end;
    
    % decode inputs
    % -------------
    %if res.plotgroups & res.plotconditions, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end;
    if res.plotgroups, res.plotgroups = 'together'; else res.plotgroups = 'apart'; end;
    if res.plotconditions , res.plotconditions  = 'together'; else res.plotconditions  = 'apart'; end;
    res.topotime  = str2num( res.topotime );
    res.timerange = str2num( res.timerange );
    res.ylim      = str2num( res.ylim );
    res.filter    = str2num( res.filter );
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( char(res.filter), char(STUDY.etc.erpparams.filter)), options = { options{:} 'filter' res.filter }; end;
    if ~strcmpi( res.plotgroups, STUDY.etc.erpparams.plotgroups), options = { options{:} 'plotgroups' res.plotgroups }; end;
    if ~strcmpi( res.plotconditions , STUDY.etc.erpparams.plotconditions ), options = { options{:} 'plotconditions'  res.plotconditions  }; end;
    if ~isequal(res.ylim     , STUDY.etc.erpparams.ylim),      options = { options{:} 'ylim' res.ylim       }; end;
    if ~isequal(res.timerange, STUDY.etc.erpparams.timerange), options = { options{:} 'timerange' res.timerange }; end;
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

