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
%  'topotime'    - [float] plot scalp map at specific time. A time range may
%                  also be provide and the erpim will be averaged over the
%                  given time range. Requires 'topofreq' below to be set.
%  'colorlimits' - [min max] color limits.
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

function [ STUDY, com ] = pop_erpimparams(STUDY, varargin);

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    icaFlag = fastif(isnan(STUDY.etc.erpimparams.topotime), 1, 0);
    
    if strcmpi(STUDY.etc.erpimparams.averagechan,'off')
        multipleChansVal   = 1; % scalp array
    else
        multipleChansVal   = 2; % average channels
    end
    
    uilist = { ...
        {'style' 'text'       'string' 'ERPimage plotting options' 'fontweight' 'bold' 'tag', 'erpim' } ...
        {'style' 'text'       'string' 'Time range in ms [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.timerange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Color limits [Low High]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.erpimparams.colorlimits) 'tag' 'colorlimits' } ...
        {} ...
        {'style' 'text'       'string' 'Multiple channels selection' 'fontweight' 'bold' 'tag', 'spec' 'fontsize', 12} ...
        {} {'style' 'popupmenu'  'string' { 'Plot channels individually' 'Average selected channels' } 'value' multipleChansVal 'tag' 'multiplechan' } };
    evalstr = 'set(findobj(gcf, ''tag'', ''erpim''), ''fontsize'', 12);';
    cbline = [0.07 1.1];
    chanline  = [ 0.07 0.8];
    geometry = { 1 [0.6 .4] [0.6 .4] 1 1 chanline };    
    
    if icaFlag
        uilist = uilist(1:end-4);
        geometry = geometry(1:end-3);
    end
    
    [~, ~, ~, res] = inputgui( 'geometry' , geometry, 'uilist', uilist, ...
                                            'title', 'Set Erpimage plotting parameters -- pop_erpimparams()', 'eval', evalstr);
    if isempty(res), return; end
    
    % decode input
    % ------------
    res.topotime    = []; %str2num( res.topotime );
    res.timerange   = str2num( res.timerange );
    res.colorlimits = str2num( res.colorlimits );
    
    % build command call
    % ------------------
    options = {};
    if ~isequal(res.topotime  ,  STUDY.etc.erpimparams.topotime),    options = { options{:} 'topotime'    res.topotime    }; end
    if ~isequal(res.timerange ,  STUDY.etc.erpimparams.timerange),   options = { options{:} 'timerange'   res.timerange   }; end
    if ~isequal(res.colorlimits, STUDY.etc.erpimparams.colorlimits), options = { options{:} 'colorlimits' res.colorlimits }; end
    if isfield(res, 'multiplechan') 
        if res.multiplechan == 1 res.multiplechan = 'off'; else res.multiplechan = 'on'; end
        if ~isequal(res.multiplechan, STUDY.etc.erpimparams.averagechan), options = { options{:} 'averagechan' res.multiplechan }; end
    end
    
    % execute options
    % ---------------
    if ~isempty(options)
        STUDY = pop_erpimparams(STUDY, options{:});
        com = sprintf('STUDY = pop_erpimparams(STUDY, %s);', vararg2str( options ));
    end
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        allfields = fieldnames(STUDY.etc.erpimparams);
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, allfields, 'exact'))
                STUDY.etc.erpimparams = setfield(STUDY.etc.erpimparams, varargin{index}, varargin{index+1});
            else
                if ~isequal(varargin{index}, 'datatype') && ~isequal(varargin{index}, 'channels')
                    inderpimopt = strmatch(varargin{index}, STUDY.etc.erpimparams.erpimageopt(1:2:end), 'exact');
                    if ~isempty(inderpimopt)
                        STUDY.etc.erpimparams.erpimageopt{2*inderpimopt} = varargin{index+1};
                    else
                        STUDY.etc.erpimparams.erpimageopt = { STUDY.etc.erpimparams.erpimageopt{:} varargin{index}, varargin{index+1} };
                    end
                end
            end
        end
    end
end

% scan clusters and channels to remove erpimdata info if timerange etc. have changed
% ---------------------------------------------------------------------------------
if ~isequal(STUDY.etc.erpimparams.timerange, TMPSTUDY.etc.erpimparams.timerange)
    rmfields = { 'erpimdata' 'erpimtimes' 'erpimtrials' 'erpimevents' };
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
    if ~isfield(STUDY.etc, 'erpimparams'), STUDY.etc.erpimparams = []; end
    if ~isfield(STUDY.etc.erpimparams, 'erpimageopt'),  STUDY.etc.erpimparams.erpimageopt = {}; end
    if ~isfield(STUDY.etc.erpimparams, 'sorttype'   ),  STUDY.etc.erpimparams.sorttype    = ''; end
    if ~isfield(STUDY.etc.erpimparams, 'sortwin'    ),  STUDY.etc.erpimparams.sortwin     = []; end
    if ~isfield(STUDY.etc.erpimparams, 'sortfield'  ),  STUDY.etc.erpimparams.sortfield   = 'latency'; end

    if ~isfield(STUDY.etc.erpimparams, 'rmcomps'    ),  STUDY.etc.erpimparams.rmcomps     = []; end
    if ~isfield(STUDY.etc.erpimparams, 'interp'     ),  STUDY.etc.erpimparams.interp      = []; end
    
    if ~isfield(STUDY.etc.erpimparams, 'timerange'  ),  STUDY.etc.erpimparams.timerange   = []; end
    if ~isfield(STUDY.etc.erpimparams, 'topotime'   ),  STUDY.etc.erpimparams.topotime    = []; end
    if ~isfield(STUDY.etc.erpimparams, 'colorlimits'),  STUDY.etc.erpimparams.colorlimits = []; end

    if ~isfield(STUDY.etc.erpimparams, 'concatenate'),  STUDY.etc.erpimparams.concatenate = 'off'; end
    if ~isfield(STUDY.etc.erpimparams, 'nlines'),       STUDY.etc.erpimparams.nlines      = 20; end
    if ~isfield(STUDY.etc.erpimparams, 'smoothing'),    STUDY.etc.erpimparams.smoothing   = 10; end
    
    if ~isfield(STUDY.etc.erpimparams, 'averagemode' ), STUDY.etc.erpimparams.averagemode  = 'ave'; end
    if ~isfield(STUDY.etc.erpimparams, 'averagechan' ), STUDY.etc.erpimparams.averagechan  = 'off'; end
