% std_erpplot() - plotting and statistics option for ERPs.
%
% Usage:    
%   >> STUDY = std_erpplot(STUDY, 'key', 'val');   
%
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG
%
% Optional inputs:
% To be documented...
%
% See also: std_erpplot()
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2006-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [ STUDY, com ] = pop_erpparams(STUDY, varargin);

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    enablegroup = fastif(length(STUDY.group)>1, 'on', 'off');
    enablecond  = fastif(length(STUDY.condition)>1, 'on', 'off');
    threshstr   = fastif(isnan(STUDY.etc.erpparams.threshold),'', num2str(STUDY.etc.erpparams.threshold));
    plotcond    = fastif(strcmpi(STUDY.etc.erpparams.plotcond, 'together'), 1, 0);
    plotgroup   = fastif(strcmpi(STUDY.etc.erpparams.plotgroup,'together'), 1, 0);
    statval     = fastif(strcmpi(STUDY.etc.erpparams.statistics,'param'), 1, 2);
    statcond    = fastif(strcmpi(STUDY.etc.erpparams.statcond, 'on'), 1, 0);
    statgroup   = fastif(strcmpi(STUDY.etc.erpparams.statgroup,'on'), 1, 0);
    
    uilist = { ...
        {'style' 'text'       'string' 'Time range (ms)'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.timerange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Plot limit (uV)'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpparams.ylim) 'tag' 'ylim' } ...
        {} {'style' 'checkbox'   'string' '' 'value' plotcond 'enable' enablecond  'tag' 'plotcond' } ...
        {'style' 'text'       'string' 'Plot conditions on the same panel' 'enable' enablecond } ...
        {} {'style' 'checkbox'   'string' '' 'value' plotgroup 'enable' enablegroup 'tag' 'plotgroup' } ...
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
                                   'helpcom', 'pophelp(''std_erpparams'')', ...
                                   'title', 'Set parameters for plotting ERPs -- pop_erpparams()');

    if isempty(res), return; end;
    
    % decode inputs
    % -------------
    if res.plotgroup & res.plotcond, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end;
    if res.statgroup, res.statgroup = 'on'; else res.statgroup = 'off'; end;
    if res.statcond , res.statcond  = 'on'; else res.statcond  = 'off'; end;
    if res.plotgroup, res.plotgroup = 'together'; else res.plotgroup = 'appart'; end;
    if res.plotcond , res.plotcond  = 'together'; else res.plotcond  = 'appart'; end;
    res.timerange = str2num( res.timerange );
    res.ylim      = str2num( res.ylim );
    res.threshold = str2num( res.threshold );
    if isempty(res.threshold),res.threshold = NaN; end;
    if res.statistics == 1, res.statistics  = 'param'; 
    else                    res.statistics  = 'perm'; 
    end;
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( res.plotgroup, STUDY.etc.erpparams.plotgroup), options = { options{:} 'plotgroup' res.plotgroup }; end;
    if ~strcmpi( res.plotcond , STUDY.etc.erpparams.plotcond ), options = { options{:} 'plotcond'  res.plotcond  }; end;
    if ~strcmpi( res.statgroup, STUDY.etc.erpparams.statgroup), options = { options{:} 'statgroup' res.statgroup }; end;
    if ~strcmpi( res.statcond , STUDY.etc.erpparams.statcond ), options = { options{:} 'statcond'  res.statcond  }; end;
    if ~strcmpi( res.statistics, STUDY.etc.erpparams.statistics ), options = { options{:} 'statistics' res.statistics }; end;
    if ~isequal(res.ylim     , STUDY.etc.erpparams.ylim),      options = { options{:} 'ylim' res.ylim      }; end;
    if ~isequal(res.timerange, STUDY.etc.erpparams.timerange), options = { options{:} 'timerange' res.timerange }; end;
    if isnan(res.threshold) & ~isnan(STUDY.etc.erpparams.threshold) | ...
            ~isnan(res.threshold) & isnan(STUDY.etc.erpparams.threshold) | ...
                ~isnan(res.threshold) & res.threshold ~= STUDY.etc.erpparams.threshold
                options = { options{:} 'threshold' res.threshold }; 
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
            STUDY.etc.erpparams = setfield(STUDY.etc.erpparams, varargin{index}, varargin{index+1});
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
    if ~isfield(STUDY.etc.erpparams, 'timerange'),  STUDY.etc.erpparams.timerange = []; end;
    if ~isfield(STUDY.etc.erpparams, 'ylim'     ),  STUDY.etc.erpparams.ylim      = []; end;
    if ~isfield(STUDY.etc.erpparams, 'statistics'), STUDY.etc.erpparams.statistics = 'param'; end;
    if ~isfield(STUDY.etc.erpparams, 'statgroup'),  STUDY.etc.erpparams.statgroup = 'off'; end;
    if ~isfield(STUDY.etc.erpparams, 'statcond' ),  STUDY.etc.erpparams.statcond  = 'off'; end;
    if ~isfield(STUDY.etc.erpparams, 'threshold' ), STUDY.etc.erpparams.threshold = NaN; end;
    if ~isfield(STUDY.etc.erpparams, 'plotgroup') , STUDY.etc.erpparams.plotgroup = 'appart'; end;
    if ~isfield(STUDY.etc.erpparams, 'plotcond') ,  STUDY.etc.erpparams.plotcond  = 'appart'; end;

