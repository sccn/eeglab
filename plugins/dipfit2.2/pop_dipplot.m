% pop_dipplot() - plot dipoles.
%
% Usage:
%   >> pop_dipplot( EEG ); % pop up interactive window
%   >> pop_dipplot( EEG, type, comps, 'key1', 'val1', 'key2', 'val2', ...);
%
% Graphic interface:
%   "Components" - [edit box] enter component number to plot. By
%                all the localized components are plotted. Command
%                line equivalent: components.
%   "Background image" - [list box] use BESA background image or average
%                MRI image. Dipplot command line equivalent: 'image'.
%   "Summary mode" - [Checkbox] when checked, plot the 3 views of the
%                head model and dipole locations.
%   "Normalized dipole length" - [Checkbox] normalize the length of
%               all dipoles. Dipplot command line equivalent: 'normlen'.
%   "Additionnal dipfit() options" - [checkbox] enter additionnal 
%               sequence of 'key', 'val' argument in this edit box.
%
% Inputs:
%   EEG   - Input dataset
%   type  - ['DIPFIT'|'BESA'] use either 'DIPFIT' dipoles or
%           'BESA' dipoles.
%   comps - [integer array] plot component indices. If empty
%           all the localized components are plotted.
%
% Optional inputs:
%   Same as dipplot().
%
% Author: Arnaud Delorme, CNL / Salk Institute, 26 Feb 2003
%
% See also: dipplot()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.4  2003/03/13 19:51:07  arno
% updated for release
%
% Revision 1.3  2003/03/11 23:35:27  arno
% typo
%
% Revision 1.2  2003/03/11 23:34:07  arno
% adding doc and options
%
% Revision 1.1  2003/02/26 17:07:14  arno
% Initial revision
%

function [com] = pop_dipplot( EEG, typedip, comps, varargin);

com ='';
if nargin < 1
    help pop_dipplot;
    return;
end;

if nargin < 3
	% popup window parameters
	% -----------------------
    varstr = '';;
    if isfield(EEG, 'dipfit'), varstr = 'DIPFIT|'; typedip = 'dipfit'; end;
    if isfield(EEG, 'sources'), varstr = [varstr 'BESA|']; typedip = 'besa';end;
    if isempty(varstr), error('No dipole information in dataset'); end;
    
    geometry = { [2 1] [2 1] [2 1] [2.05 0.23 .75] [2.05 0.23 .75] [2 1] };
    uilist = { { 'style' 'text' 'string' 'Components indices ([]=all)' } ...
               { 'style' 'edit' 'string' '' } ...
               { 'style' 'text' 'string' 'Use dipoles from (scroll, then click to select)' } ...
               { 'style' 'listbox' 'string' varstr(1:end-1) 'value' 1 } ...
               { 'style' 'text' 'string' 'Background image (click to select)' } ...
               { 'style' 'listbox' 'string' 'BESA Head|average MRI' } ...
               { 'style' 'text' 'string' 'Sumary mode' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Normalized dipole length' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Additionnal dipfit() options' } ...
               { 'style' 'edit' 'string' '' } };
     
    if length(varstr) < 7,
        uilist(3:4) = []; geometry(2) = [];
    end;
	result = inputgui( geometry, uilist, 'pophelp(''pop_dipplot'')', 'Plot dipoles - pop_dipplot');
	if length(result) == 0 return; end;

	% decode parameters
	% -----------------
    options = {};
    if ~isempty(result{1}), comps = eval( [ '[' result{1} ']' ] ); else comps = []; end;
    if length(varstr) >= 7,
        if result{2} == 1, typedip = 'DIPFIT';
        else               typedip = 'BESA';
        end;
        ind = 3;
    else
        ind = 2;
    end;
    options = { options{:} 'image' fastif(result{ind} == 2, 'mri', 'besa') };
    if result{ind+1} == 1, options = { options{:} 'summary' 'on' }; end;
    if result{ind+2} == 1, options = { options{:} 'normlen' 'on' }; end;
    if ~isempty( result{ind+3} ), options = { options{:} eval( [ '{' result{ind+3} '}' ] ) }; end;
else 
    options = varargin;
end;

if strcmpi(typedip, 'besa')
    if ~isfield(EEG, 'sources'), error('No BESA dipole information in dataset');end;
    if ~isempty(comps)
        [tmp1 int] = intersect(cell2mat({EEG.sources.component}), comps);
        if isempty(int), error ('Localization not found for selected components'); end;
        dipplot(EEG.sources(int), options{:});
    else
        dipplot(EEG.sources, options{:});
    end;      
else 
    if ~isfield(EEG, 'dipfit'), error('No DIPFIT dipole information in dataset');end;
    if ~isempty(comps)
        dipplot(EEG.dipfit.model(comps), options{:});
    else
        dipplot(EEG.dipfit.model, options{:});
    end;
end;
    
com = sprintf('pop_dipplot( %s, ''%s'', %s);', inputname(1), typedip, vararg2str(options));
return;
