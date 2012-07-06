% pop_erpimparams() - Set plotting and statistics parameters for 
%                    computing and plotting STUDY mean ERPimages and measure 
%                    statistics. Settings are stored within the STUDY 
%                    structure (STUDY.etc.erpimparams) which is used
%                    whenever plotting is performed by the function
%                    std_erpimage().
% Usage:    
%   >> STUDY = pop_erpimparams(STUDY, 'key', 'val', ...);   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% Erpimage plotting options:
%  'timerange'   - [min max] erpim/ITC plotting latency range in ms. 
%                  {default: the whole output latency range}.
%  'trialrange'  - [min max] erpim/ITC plotting frequency range in ms. 
%                  {default: the whole output frequency range}
%  'topotime'    - [float] plot scalp map at specific time. A time range may
%                  also be provide and the erpim will be averaged over the
%                  given time range. Requires 'topofreq' below to be set.
%  'topotrial'  - [float] plot scalp map at specific trial in ERPimage. As 
%                  above a trial range may also be provided.
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

function [ STUDY, com ] = pop_erpimparams(STUDY, varargin);

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    vis = fastif(isnan(STUDY.etc.erpimparams.topotime), 'off', 'on');
    uilist = { ...
        {'style' 'text'       'string' 'ERPimage plotting options' 'fontweight' 'bold' 'tag', 'erpim' } ...
        {'style' 'text'       'string' 'Time range in ms [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.timerange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Plot scalp map at time [ms]' 'visible' vis} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.topotime) 'tag' 'topotime' 'visible' vis } ...
        {'style' 'text'       'string' 'Trial range [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.trialrange) 'tag' 'trialrange' } ...
        {'style' 'text'       'string' 'Plot scalp map at trial(s)' 'visible' vis} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.topotrial) 'tag' 'topotrial' 'visible' vis } ...
        {'style' 'text'       'string' 'Color limits [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.colorlimits) 'tag' 'colorlimits' } ...
        {'style' 'text'       'string' '' } ...
        {'style' 'text'       'string' '' } };
    evalstr = 'set(findobj(gcf, ''tag'', ''erpim''), ''fontsize'', 12);';
    cbline = [0.07 1.1];
    otherline = [ 0.6 .4 0.6 .4];
    geometry = { 1 otherline otherline otherline };
    enablecond  = fastif(length(STUDY.design(STUDY.currentdesign).variable(1).value)>1, 'on', 'off');
    enablegroup = fastif(length(STUDY.design(STUDY.currentdesign).variable(2).value)>1, 'on', 'off');
    
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, ...
                                            'title', 'Set Erpimage plotting parameters -- pop_erpimparams()', 'eval', evalstr);
    if isempty(res), return; end;
    
    % decode input
    % ------------
    res.topotime    = str2num( res.topotime );
    res.topotrial   = str2num( res.topotrial );
    res.timerange   = str2num( res.timerange );
    res.trialrange  = str2num( res.trialrange );
    res.colorlimits = str2num( res.colorlimits );
    
    % build command call
    % ------------------
    options = {};
    if ~isequal(res.topotime  ,  STUDY.etc.erpimparams.topotime),    options = { options{:} 'topotime'    res.topotime    }; end;
    if ~isequal(res.topotrial,   STUDY.etc.erpimparams.topotrial),   options = { options{:} 'topotrial'   res.topotrial   }; end;
    if ~isequal(res.timerange ,  STUDY.etc.erpimparams.timerange),   options = { options{:} 'timerange'   res.timerange   }; end;
    if ~isequal(res.trialrange,  STUDY.etc.erpimparams.trialrange),  options = { options{:} 'trialrange'  res.trialrange  }; end;
    if ~isequal(res.colorlimits, STUDY.etc.erpimparams.colorlimits), options = { options{:} 'colorlimits' res.colorlimits }; end;
    if ~isempty(options)
        STUDY = pop_erpimparams(STUDY, options{:});
        com = sprintf('STUDY = pop_erpimparams(STUDY, %s);', vararg2str( options ));
    end;
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        allfields = fieldnames(STUDY.etc.erpimparams);
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, allfields, 'exact'))
                STUDY.etc.erpimparams = setfield(STUDY.etc.erpimparams, varargin{index}, varargin{index+1});
            else
                inderpimopt = strmatch(varargin{index}, STUDY.etc.erpimparams.erpimageopt(1:2:end), 'exact');
                if ~isempty(inderpimopt)
                    STUDY.etc.erpimparams.erpimageopt{inderpimopt+1} = varargin{index+1};
                else
                    STUDY.etc.erpimparams.erpimageopt = { STUDY.etc.erpimparams.erpimageopt{:} varargin{index}, varargin{index+1} };
                end;
            end;
        end;
    end;
end;

% scan clusters and channels to remove erpimdata info if timerange etc. have changed
% ---------------------------------------------------------------------------------
if ~isequal(STUDY.etc.erpimparams.timerange, TMPSTUDY.etc.erpimparams.timerange) | ... 
    ~isequal(STUDY.etc.erpimparams.trialrange, TMPSTUDY.etc.erpimparams.trialrange)
    if isfield(STUDY.cluster, 'erpimdata')
        for index = 1:length(STUDY.cluster)
            STUDY.cluster(index).erpimdata   = [];
            STUDY.cluster(index).erpimtimes  = [];
            STUDY.cluster(index).erpimtrials = [];
            STUDY.cluster(index).erpimevents = [];
        end;
    end;
    if isfield(STUDY.changrp, 'erpimdata')
        for index = 1:length(STUDY.changrp)
            STUDY.changrp(index).erpimdata   = [];
            STUDY.changrp(index).erpimtimes  = [];
            STUDY.changrp(index).erpimtrials = [];
            STUDY.changrp(index).erpimevents = [];
        end;
    end;
end;

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'erpimparams'), STUDY.etc.erpimparams = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'erpimageopt'),  STUDY.etc.erpimparams.erpimageopt = {}; end;
    if ~isfield(STUDY.etc.erpimparams, 'sorttype'   ),  STUDY.etc.erpimparams.sorttype    = ''; end;
    if ~isfield(STUDY.etc.erpimparams, 'sortwin'    ),  STUDY.etc.erpimparams.sortwin     = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'sortfield'  ),  STUDY.etc.erpimparams.sortfield   = 'latency'; end;

    if ~isfield(STUDY.etc.erpimparams, 'rmcomps'    ),  STUDY.etc.erpimparams.rmcomps     = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'interp'     ),  STUDY.etc.erpimparams.interp      = []; end;
    
    if ~isfield(STUDY.etc.erpimparams, 'timerange'  ),  STUDY.etc.erpimparams.timerange   = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'trialrange' ),  STUDY.etc.erpimparams.trialrange  = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'topotime'   ),  STUDY.etc.erpimparams.topotime    = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'topotrial' ),   STUDY.etc.erpimparams.topotrial   = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'colorlimits'),  STUDY.etc.erpimparams.colorlimits = []; end;

    if ~isfield(STUDY.etc.erpimparams, 'concatenate'),  STUDY.etc.erpimparams.concatenate = 'off'; end;
    if ~isfield(STUDY.etc.erpimparams, 'nlines'),       STUDY.etc.erpimparams.nlines      = 20; end;
    if ~isfield(STUDY.etc.erpimparams, 'smoothing'),    STUDY.etc.erpimparams.smoothing   = 10; end;
