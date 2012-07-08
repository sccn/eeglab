% pop_eeglab2loreta() - export EEGLAB data and channel locations to LORETA
%
% Usage:
%   >> pop_eeglab2loreta( EEG, 'key', 'val' );
%
% Inputs:
%   EEG            - EEGLAB dataset structure
%
% Optional inputs:
%   'excludechans' - indices of channels to exclude
%   'labelonly'    - only export channel labels (which position can be
%                    be then looked up in LORETA)
%   'compnum'      - indices of components.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2005
%
% See also: eeglab2loreta()

% Copyright (C) 2005 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function command = pop_eeglab2loreta(EEG, varargin); 

    command = '';
    if nargin < 1
        help pop_eeglab2loreta;
        return;
    end;    
    
    if nargin > 2
        eeglab2loreta(EEG.chanlocs, EEG.icawinv, varargin{:});
        return;
    end;
    
    commandload = [ 'filepath = uigetdir(pwd, ''Select a folder'');' ...
                    'if filepath ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''folder''), ''string'', filepath);' ...
                    'end;' ...
                    'clear filepath;' ];
    cb_plotcomps = [ 'compinds = str2num(get(findobj(gcbf, ''tag'', ''compind''), ''string''));' ...
                     'pop_topoplot(EEG, 0, compinds);' ];
    cb_erp = 'set(findobj(gcbf, ''tag'', ''compind''), ''enable'', fastif(isempty(get(gcbo, ''string'')), ''on'', ''off''));';
                 
    listui = { { 'style' 'text'       'string' 'Ouput folder' } ...
               { 'style' 'edit'       'string' '' 'tag' 'folder' } ...
               { 'style' 'pushbutton' 'string' 'Browse' 'callback' commandload } ...
               { 'style' 'text'       'string' 'Export channel labels only' } ...
               { 'style' 'checkbox'   'string' '' } ...
               { 'style' 'text'       'string' '(check=yes)' } ...
               { 'style' 'text'       'string' 'Omit channel indices' } ...
               { 'style' 'edit'       'string' '' } ...
               { } ...
               { 'style' 'text'       'string' 'Export ICA component indices' } ...
               { 'style' 'edit'       'string' [ '1:' num2str(size(EEG.icawinv,2)) ] 'tag', 'compind' } ...
               { 'style' 'pushbutton' 'string' 'Plot topo.' 'callback' cb_plotcomps } ...
               { 'style' 'text'       'string' 'or export ERP time range [min max] (ms)' } ...
               { 'style' 'edit'       'string' '' 'tag', 'erp' 'callback' cb_erp } ...
               {  } ...
               };
    geom = { [1 0.5 0.5] [1 0.5 0.5] [1 0.5 0.5] [1 0.5 0.5] [1 0.5 0.5] };
    results = inputgui('geometry', geom, 'uilist', listui, 'helpcom', 'pophelp(''pop_eeglab2loreta'')', ...
        'title', 'Export EEG info to LORETA');
    if isempty(results), return; end;
        
    % decode inputs
    % -------------
    folderout = results{1};
    options = {};
    if results{2}
        options = { 'labelonly' 'on' };
    end;
    if ~isempty(results{3})
        options = { options{:} 'excludechan' eval( [ '[' results{3} ']' ] ) };
    end;
    
    % export erp
    % ----------
    if ~isempty(results{5})
        timerange = eval( [ '[' results{5} ']' ] );
        [tmp minind1] = min(abs(EEG.times-timerange(1)));
        [tmp minind2] = min(abs(EEG.times-timerange(2)));
        eeglab2loreta(EEG.chanlocs, mean(EEG.data(:, minind1:minind2, :), 3), 'exporterp', 'on', options{:});
        return;
    end;
        
    % export comps
    % ------------
    if ~isempty(results{4})
        compnums = eval( [ '[' results{4} ']' ] );
    end;
    eeglab2loreta(EEG.chanlocs, EEG.icawinv, 'compnum', compnums, options{:});
