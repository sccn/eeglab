% pop_dipplot() - plot dipoles.
%
% Usage:
%   >> pop_dipplot( EEG ); % pop up interactive window
%   >> pop_dipplot( EEG, comps, 'key1', 'val1', 'key2', 'val2', ...);
%
% Graphic interface:
%   "Components" - [edit box] enter component number to plot. By
%                all the localized components are plotted. Command
%                line equivalent: components.
%   "Background image" - [edit box] MRI background image. This image 
%               has to be normalized to the MNI brain using SPM2 for
%               instance. Dipplot() command line equivalent: 'image'.
%   "Summary mode" - [Checkbox] when checked, plot the 3 views of the
%                head model and dipole locations. Dipplot() equivalent 
%               is 'summary' and 'num'.
%   "Plot edges" - [Checkbox] plot edges at the intersection between
%               MRI slices. Diplot() equivalent is 'drawedges'.
%   "Plot closest MRI slide" - [Checkbox] plot closest MRI slice to
%               dipoles although not using the 'tight' view mode.
%               Dipplot() equivalent is 'cornermri' and 'axistight'.
%   "Plot dipole's 2-D projections" - [Checkbox] plot a dimed dipole
%               projection on each 2-D MRI slice. Dipplot() equivalent 
%               is 'projimg'.
%   "Plot projection lines" - [Checkbox] plot lines originating from
%               dipoles and perpendicular to each 2-D MRI slice. 
%               Dipplot() equivalent is 'projline'.
%   "Make all dipole point out" - [Checkbox] make all dipole point 
%               toward outside the brain. Dipplot() equivalent is 
%               'pointout'.
%   "Normalized dipole length" - [Checkbox] normalize the length of
%               all dipoles. Dipplot() command line equivalent: 'normlen'.
%   "Additionnal dipfit() options" - [checkbox] enter additionnal 
%               sequence of 'key', 'val' argument in this edit box.
%
% Inputs:
%   EEG   - Input dataset
%   comps - [integer array] plot component indices. If empty
%           all the localized components are plotted.
%
% Optional inputs:
%   'key','val' - same as dipplot()
%
% Author: Arnaud Delorme, CNL / Salk Institute, 26 Feb 2003-
%
% See also: dipplot()

%   "Use dipoles from" - [list box] use dipoles from BESA or from the
%                DIPFIT toolbox. Command line equivalent: type.

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [com] = pop_dipplot( EEG, comps, varargin);

com ='';
if nargin < 1
    help pop_dipplot;
    return;
end;

% check input structure
% ---------------------
if ~isfield(EEG, 'dipfit') & ~isfield(EEG, 'sources')
    if ~isfield(EEG.dipfit.hdmfile) & ~isfield(EEG, 'sources')
        error('No dipole information in dataset'); 
    end;
    error('No dipole information in dataset'); 
end;
if ~isfield(EEG.dipfit, 'model')
    error('No dipole information in dataset'); 
end;
    
typedip = 'nonbesa';
if nargin < 2
	% popup window parameters
	% -----------------------
    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''mrifile''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    
    geometry = { [2 1] [2 1] [0.8 0.3 1.5] [2.05 0.26 .75] [2.05 0.26 .75] [2.05 0.26 .75] ...
                 [2.05 0.26 .75] [2.05 0.26 .75] [2.05 0.26 .75] [2.05 0.26 .75] [2 1] };
    uilist = { { 'style' 'text' 'string' 'Components indices ([]=all avaliable)' } ...
               { 'style' 'edit' 'string' '' } ...
               { 'style' 'text' 'string' 'Plot dipoles within RV (%) range ([min max])' } ...
               { 'style' 'edit' 'string' '' } ...
               { 'style' 'text' 'string' 'Background image' } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' commandload } ...
               { 'style' 'edit' 'string' EEG.dipfit.mrifile 'tag' 'mrifile' } ...
               { 'style' 'text' 'string' 'Plot summary mode' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Plot edges' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Plot closest MRI slide' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Plot dipole''s 2-D projections' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Plot projection lines' } ...
               { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
               { 'style' 'text' 'string' 'Make all dipoles point out' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Normalized dipole length' } ...
               { 'style' 'checkbox' 'string' '' 'value' 1 } {} ...
               { 'style' 'text' 'string' 'Additionnal dipplot() options' } ...
               { 'style' 'edit' 'string' '' } };
     
	result = inputgui( geometry, uilist, 'pophelp(''pop_dipplot'')', 'Plot dipoles - pop_dipplot');
	if length(result) == 0 return; end;

	% decode parameters
	% -----------------
    options = {};
    if ~isempty(result{1}), comps = eval( [ '[' result{1} ']' ] ); else comps = []; end;
    if ~isempty(result{2}), options = { options{:} 'rvrange' eval(  [ '[' result{2} ']' ] ) }; end;
    options = { options{:} 'mri' result{3} };
    if result{4} == 1, options = { options{:} 'summary'   'on' 'num' 'on' }; end;
    if result{5} == 1, options = { options{:} 'drawedges' 'on' }; end;
    if result{6} == 1, options = { options{:} 'cornermri' 'on' 'axistight' 'on' }; end;
    if result{7} == 1, options = { options{:} 'projimg'   'on' }; end;
    if result{8} == 1, options = { options{:} 'projlines' 'on' }; end;
    if result{9} == 1, options = { options{:} 'pointout'  'on' }; end; 
    if result{10} == 1, options = { options{:} 'normlen'   'on' }; end;
    if ~isempty( result{11} ), tmpopt = eval( [ '{' result{11} '}' ] ); options = { options{:} tmpopt{:} }; end;
else 
    if isstr(comps)
        typedip = comps;
        options = varargin(2:end);
        comps = varargin{1};
    else
        options = varargin;
    end;
end;

if strcmpi(typedip, 'besa')
    if ~isfield(EEG, 'sources'), error('No BESA dipole information in dataset');end;
    if ~isempty(comps)
        [tmp1 int] = intersect( [ EEG.sources.component ], comps);
        if isempty(int), error ('Localization not found for selected components'); end;
        dipplot(EEG.sources(int), 'sphere', 1, options{:});
    else
        dipplot(EEG.sources, options{:});
    end;      
else 
    if ~isfield(EEG, 'dipfit'), error('No DIPFIT dipole information in dataset');end;

    % components to plot
    % ------------------
    if ~isempty(comps)
        if ~isfield(EEG.dipfit.model, 'component')
            for index = double(comps(:)')
                EEG.dipfit.model(index).component = index;
            end;
        end;
    else
        % find localized dipoles
        comps = [];
        for index2 = 1:length(EEG.dipfit.model)
            if ~isempty(EEG.dipfit.model(index2).posxyz) ~= 0
                comps = [ comps index2 ];
                EEG.dipfit.model(index2).component = index2;
            end;
        end;        
    end;
    
    % plotting
    % --------
    tmpoptions = { options{:} 'coordformat', EEG.dipfit.coordformat };
    if strcmpi(EEG.dipfit.coordformat, 'spherical')
        dipplot(EEG.dipfit.model(comps), tmpoptions{:});
    elseif strcmpi(EEG.dipfit.coordformat, 'CTF')
        dipplot(EEG.dipfit.model(comps), tmpoptions{:});
    else
        dipplot(EEG.dipfit.model(comps), 'meshdata', EEG.dipfit.hdmfile, tmpoptions{:});
    end;
end;
    
if nargin < 3
    com = sprintf('pop_dipplot( EEG, %s);', vararg2str({ comps options{:}}));
end;
return;
