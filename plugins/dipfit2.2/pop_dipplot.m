% pop_dipplot() - plot dipoles.
%
% Usage:
%   >> pop_dipplot( EEG ); % pop up interactive window
%   >> pop_dipplot( EEG, type, 'key1', 'val1', 'key2', 'val2', ...);
%
% Graphic interface:
%   
% Inputs:
%   EEG        - Input dataset
%   type       - ['DIPFIT'|'BESA'] use either 'DIPFIT' dipoles or
%                'BESA' dipoles.
%
% Optional inputs:
%   Same as dipplot().
%
% Author: Arnaud Delorme, CNL / Salk Institute, 26 Feb 2003
%
% See also: dipplot()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

function [com] = pop_dipplot( EEG, typedip, varargin);

com ='';
if nargin < 1
    help pop_dipplot;
    return;
end;

if nargin < 2
	% popup window parameters
	% -----------------------
    if isfield(EEG, 'dipfit'), defaulttype = 1;
    elseif isfield(EEG, 'sources'), defaulttype = 2;
    else error('No dipole information in dataset'); 
    end;
    geometry = { [2 1] [2 1] [2.05 0.23 .75] [2 1] };
    uilist = { { 'style' 'text' 'string' 'Use dipoles from (scroll, then click to select)' } ...
               { 'style' 'listbox' 'string' 'DIPFIT|BESA' 'value' defaulttype } ...
               { 'style' 'text' 'string' 'Background image' } ...
               { 'style' 'listbox' 'string' 'average MRI|BESA Head' } ...
               { 'style' 'text' 'string' 'Sumary mode' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Additionnal dipfit() options' } ...
               { 'style' 'edit' 'string' '' } };
               
	result = inputgui( geometry, uilist, 'pophelp(''pop_dipplot'')', 'Plot dipoles - pop_dipplot');
	if length(result) == 0 return; end;

	% decode parameters
	% -----------------
    options = {};
    if result{1} == 1, typedip = 'DIPFIT';
    else               typedip = 'BESA';
    end;
    options = { options{:} 'image' fastif(result{2} == 1, 'mri', 'besa') };
    if result{1} == 1, options = { options{:} 'summary' 'on' }; end;
    if ~isempty( result{4} ), options = { options{:} eval( [ '{' result{4} '}' ] ) }; end;
else 
    options = varargin;
end;

if strcmpi(typedip, 'besa')
    if ~isfield(EEG, 'sources'), error('No BESA dipole information in dataset');end;
    dipplot(EEG.sources, options{:});
else 
    if ~isfield(EEG, 'dipfit'), error('No DIPFIT dipole information in dataset');end;
    dipplot(EEG.sources, options{:});
end;
    
com = sprintf('pop_dipplot( %s, ''%s'', %s);', inputname(1), typedip, vararg2str(options));
return;
