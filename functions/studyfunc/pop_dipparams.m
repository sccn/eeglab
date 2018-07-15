% pop_dipparams() - Set plotting parameters for dipoles.
%
% Usage:    
%   >> STUDY = pop_dipparams(STUDY, 'key', 'val');   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% Optional inputs:
%   'axistight' - ['on'|'off'] Plot closest MRI slide. Default is 'off'.
%   'projimg'   - ['on'|'off'] lot dipoles projections on each axix. Default is 'off'.
%   'projlines' - ['on'|'off'] Plot projection lines. Default is 'off'.
%   'centrline' - ['on'|'off'] Plot centroid's projection line. 
%                 This option overwrite 'projlines'. Default is 'on'.
%
% See also: std_dipplot()
%
% Authors: Arnaud Delorme

% Copyright (C) Arnaud Delorme, 20013
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

function [ STUDY, com ] = pop_dipparams(STUDY, varargin);

STUDY = default_params(STUDY);
com = '';
if isempty(varargin)
    
    val_axistight = fastif(strcmpi(STUDY.etc.dipparams.axistight,'on'), 1, 0);
    val_projimg   = fastif(strcmpi(STUDY.etc.dipparams.projimg,'on'), 1, 0);
    val_projlines = fastif(strcmpi(STUDY.etc.dipparams.projlines,'on'), 1, 0);
    val_density   = fastif(strcmpi(STUDY.etc.dipparams.density,'on'), 1, 0);
    val_centrline = fastif(strcmpi(STUDY.etc.dipparams.centrline,'on'), 1, 0);
   
    if val_projlines, val_centrlineenable ='off'; else, val_centrlineenable ='on'; end
    
    uilist = { ...
        {'style' 'checkbox' 'tag' 'density'   'value' val_density   } { 'style' 'text'  'string' 'Plot dipoles density'         } ...
        {'style' 'checkbox' 'tag' 'projlines' 'value' val_projlines...
        'callback' ['set(findobj(gcf,''tag'',''centrline''),''value'',get(findobj(gcf,''tag'',''projlines''),''value''));'...
         'if get(findobj(gcf,''tag'',''projlines''),''value''),set(findobj(gcf,''tag'',''centrline''),''Enable'',''off'');'...
         ' else,set(findobj(gcf,''tag'',''centrline''),''Enable'',''on''); end;'] }...
        { 'style' 'text'  'string' 'Plot projection lines'} ...
        {} {'style' 'checkbox' 'tag' 'centrline' 'Enable' val_centrlineenable 'value' val_centrline } { 'style' 'text'  'string' 'Plot centroid projection line'} ...
        {'style' 'checkbox' 'tag' 'projimg'   'value' val_projimg   } { 'style' 'text'  'string' 'Plot dipoles projections'     } ...
        {'style' 'checkbox' 'tag' 'axistight' 'value' val_axistight } { 'style' 'text'  'string' 'Plot closest MRI slide'       } ...
        };

    [out_param, userdat, tmp, res] = inputgui( 'geometry' , { [0.2 1]  [0.2 1] [0.1 0.2 1] [0.2 1] [0.2 1] }, 'uilist', uilist, 'geomvert', [1 1 1 1 1], ...
                                            'title', 'ERP plotting options -- pop_dipparams()');
    if isempty(res), return; end
    
    % decode inputs
    % -------------
    res.axistight = fastif(res.axistight, 'on', 'off');
    res.projimg   = fastif(res.projimg  , 'on', 'off');
    res.projlines = fastif(res.projlines, 'on', 'off');
    res.centrline = fastif(res.centrline, 'on', 'off');
    res.density   = fastif(res.density,   'on', 'off');
    
    % build command call
    % ------------------
    options = {};
    if ~strcmpi( res.axistight, STUDY.etc.dipparams.axistight), options = { options{:} 'axistight' res.axistight }; end
    if ~strcmpi( res.projimg,   STUDY.etc.dipparams.projimg  ), options = { options{:} 'projimg'   res.projimg   }; end
    if ~strcmpi( res.projlines, STUDY.etc.dipparams.projlines), options = { options{:} 'projlines' res.projlines }; end
    if ~strcmpi( res.centrline, STUDY.etc.dipparams.centrline), options = { options{:} 'centrline' res.centrline }; end
    if ~strcmpi( res.density  , STUDY.etc.dipparams.density),   options = { options{:} 'density'   res.density   }; end
    if ~isempty(options)
        STUDY = pop_dipparams(STUDY, options{:});
        com = sprintf('STUDY = pop_dipparams(STUDY, %s);', vararg2str( options ));
    end
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_params(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.dipparams), 'exact'))
                STUDY.etc.dipparams = setfield(STUDY.etc.dipparams, varargin{index}, varargin{index+1});
            end
        end
    end
end

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'dipparams'), STUDY.etc.dipparams = []; end
    if ~isfield(STUDY.etc.dipparams, 'axistight'),        STUDY.etc.dipparams.axistight = 'off'; end
    if ~isfield(STUDY.etc.dipparams, 'projimg'),          STUDY.etc.dipparams.projimg   = 'off'; end
    if ~isfield(STUDY.etc.dipparams, 'projlines'),        STUDY.etc.dipparams.projlines = 'off'; end
    if ~isfield(STUDY.etc.dipparams, 'density'),          STUDY.etc.dipparams.density   = 'off'; end
    if ~isfield(STUDY.etc.dipparams, 'centrline'),        STUDY.etc.dipparams.centrline = 'on';  end
    if fastif(strcmpi(STUDY.etc.dipparams.projlines,'on'), 1, 0),STUDY.etc.dipparams.centrline = 'on'; end
