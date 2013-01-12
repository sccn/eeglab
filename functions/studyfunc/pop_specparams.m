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
%   'averagechan' - ['on'|'off'] average data channels when several are
%                  selected.
%
% See also: std_specplot()
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2006-

% Copyright (C) Arnaud Delorme, CERCO
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

function [ STUDY, com ] = pop_specparams(STUDY, varargin);

STUDY = default_params(STUDY);
TMPSTUDY = STUDY;
com = '';
if isempty(varargin)
    
    enablecond  = fastif(length(STUDY.design(STUDY.currentdesign).variable(1).value)>1, 'on', 'off');
    enablegroup = fastif(length(STUDY.design(STUDY.currentdesign).variable(2).value)>1, 'on', 'off');
    plotconditions    = fastif(strcmpi(STUDY.etc.specparams.plotconditions, 'together'), 1, 0);
    plotgroups   = fastif(strcmpi(STUDY.etc.specparams.plotgroups,'together'), 1, 0);
    submean      = fastif(strcmpi(STUDY.etc.specparams.subtractsubjectmean,'on'), 1, 0);
    radio_averagechan  = fastif(strcmpi(STUDY.etc.specparams.averagechan,'on'), 1, 0);
    radio_scalptopo    = fastif(isempty(STUDY.etc.specparams.topofreq), 0, 1);
    if radio_scalptopo, radio_averagechan = 0; end;
    if radio_scalptopo+radio_averagechan == 0, radio_scalparray = 1; else radio_scalparray = 0; end;
        
    cb_radio = [ 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0);' ...
                 'set(gcbo, ''value'', 1);' ...
                 'set(findobj(gcbf, ''tag'', ''topofreq''), ''string'', '''');' ];
    cb_edit  = [ 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0);' ...
                 'set(findobj(gcbf, ''tag'', ''scalptopotext''), ''value'', 1);' ];
    
    uilist = { ...
        {'style' 'text'       'string' 'Spectrum plotting options' 'fontweight' 'bold' 'fontsize', 12} ...
        {} {'style' 'text'       'string' 'Frequency [low_Hz high_Hz]' } ...
        {'style' 'edit'       'string' num2str(STUDY.etc.specparams.freqrange) 'tag' 'freqrange' } ...
        {} {'style' 'text'       'string' 'Plot limits [low high]'} ...
        {'style' 'edit'       'string' num2str(STUDY.etc.specparams.ylim) 'tag' 'ylim' } ...
        {} {'style' 'checkbox'   'string' 'Subtract individual subject mean spectrum' 'value' submean 'tag' 'submean' } ...
        {} ...
        {'style' 'text'       'string' 'Spectrum plotting format' 'fontweight' 'bold' 'fontsize', 12} ...
        {} {'style' 'checkbox'   'string' 'Plot first variable on the same panel' 'value' plotconditions 'enable' enablecond  'tag' 'plotconditions' } ...
        {} {'style' 'checkbox'   'string' 'Plot second variable on the same panel' 'value' plotgroups 'enable' enablegroup 'tag' 'plotgroups' } ...
        {} ...
        {'style' 'text'       'string' 'Multiple channels selection' 'fontweight' 'bold' 'fontsize', 12} ...
        {} {'style' 'radio'   'string' 'Plot channels in scalp array'    'value' radio_scalparray 'tag' 'scalparray'       'userdata' 'radio' 'callback' cb_radio} { } ...
        {} {'style' 'radio'   'string'  'Plot topography at freq. (Hz)' 'value' radio_scalptopo  'tag' 'scalptopotext' 'userdata' 'radio' 'callback' cb_radio} ...
           {'style' 'edit'    'string' num2str(STUDY.etc.specparams.topofreq) 'tag' 'topofreq' 'callback' cb_edit } ...
        {} {'style' 'radio'   'string' 'Average selected channels' 'value' radio_averagechan 'tag' 'averagechan' 'userdata' 'radio' 'callback' cb_radio} { } };
    cbline = [0.07 1.1];
    otherline = [ 0.07 0.6 .3];
    chanline  = [ 0.07 0.8 0.3];
    geometry = { 1 otherline otherline cbline 1 1 cbline cbline 1 1 chanline chanline chanline };
    geomvert = [1.2 1 1 1 0.5 1.2 1 1 0.5 1.2 1 1 1 ];
    
    % component plotting
    % ------------------
    if isnan(STUDY.etc.specparams.topofreq)
        geometry(end-4:end) = []; 
        geomvert(end-4:end) = []; 
        uilist(end-10:end) = [];
    end;
    
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'geomvert', geomvert, ...
                                            'title', 'Spectrum plotting options -- pop_specparams()');
    if isempty(res), return; end;
    
    % decode inputs
    % -------------
    %if res.plotgroups & res.plotconditions, warndlg2('Both conditions and group cannot be plotted on the same panel'); return; end;
    if res.submean   , res.submean    = 'on'; else res.submean    = 'off'; end;
    if res.plotgroups, res.plotgroups = 'together'; else res.plotgroups = 'apart'; end;
    if res.plotconditions , res.plotconditions  = 'together'; else res.plotconditions  = 'apart'; end;
    if ~isfield(res, 'topofreq'), res.topofreq = STUDY.etc.specparams.topofreq;
    else res.topofreq  = str2num( res.topofreq );
    end;
    if ~isfield(res, 'averagechan'), res.averagechan = STUDY.etc.specparams.averagechan;
    elseif res.averagechan, res.averagechan = 'on'; else res.averagechan = 'off';
    end;
    res.freqrange = str2num( res.freqrange );
    res.ylim      = str2num( res.ylim );
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( res.plotgroups, STUDY.etc.specparams.plotgroups), options = { options{:} 'plotgroups' res.plotgroups }; end;
    if ~strcmpi( res.plotconditions , STUDY.etc.specparams.plotconditions ), options = { options{:} 'plotconditions'  res.plotconditions  }; end;
    if ~strcmpi( res.submean   , STUDY.etc.specparams.subtractsubjectmean ), options = { options{:} 'subtractsubjectmean'  res.submean  }; end;
    if ~isequal(res.topofreq, STUDY.etc.specparams.topofreq),       options = { options{:} 'topofreq' res.topofreq }; end;
    if ~isequal(res.ylim, STUDY.etc.specparams.ylim),               options = { options{:} 'ylim' res.ylim      }; end;
    if ~isequal(res.freqrange, STUDY.etc.specparams.freqrange),     options = { options{:} 'freqrange' res.freqrange }; end;
    if ~isequal(res.averagechan, STUDY.etc.specparams.averagechan), options = { options{:} 'averagechan' res.averagechan }; end;
    if ~isempty(options)
        STUDY = pop_specparams(STUDY, options{:});
        com = sprintf('STUDY = pop_specparams(STUDY, %s);', vararg2str( options ));
    end;
else
    
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.specparams), 'exact'))
                STUDY.etc.specparams = setfield(STUDY.etc.specparams, varargin{index}, varargin{index+1});
            end;
        end;
    end;
end;

% scan clusters and channels to remove specdata info if freqrange has changed
% ----------------------------------------------------------
if ~isequal(STUDY.etc.specparams.freqrange, TMPSTUDY.etc.specparams.freqrange) | ...
    ~isequal(STUDY.etc.specparams.subtractsubjectmean, TMPSTUDY.etc.specparams.subtractsubjectmean)
    rmfields = { 'specdata' 'specfreqs' };
    for iField = 1:length(rmfields)
        if isfield(STUDY.cluster, rmfields{iField})
            STUDY.cluster = rmfield(STUDY.cluster, rmfields{iField});
        end;
        if isfield(STUDY.changrp, rmfields{iField})
            STUDY.changrp = rmfield(STUDY.changrp, rmfields{iField});
        end;
    end;   
end;

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'specparams'), STUDY.etc.specparams = []; end;
    if ~isfield(STUDY.etc.specparams, 'topofreq'),   STUDY.etc.specparams.topofreq = []; end;
    if ~isfield(STUDY.etc.specparams, 'freqrange'),  STUDY.etc.specparams.freqrange = []; end;
    if ~isfield(STUDY.etc.specparams, 'ylim'     ),  STUDY.etc.specparams.ylim      = []; end;
    if ~isfield(STUDY.etc.specparams, 'subtractsubjectmean' ), STUDY.etc.specparams.subtractsubjectmean  = 'off'; end;
    if ~isfield(STUDY.etc.specparams, 'plotgroups'), STUDY.etc.specparams.plotgroups = 'apart'; end;
    if ~isfield(STUDY.etc.specparams, 'plotconditions'),  STUDY.etc.specparams.plotconditions  = 'apart'; end;
    if ~isfield(STUDY.etc.specparams, 'averagechan') ,    STUDY.etc.specparams.averagechan  = 'off'; end;
