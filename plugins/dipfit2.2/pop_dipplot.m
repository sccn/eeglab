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
% Revision 1.20  2004/03/24 01:26:20  arno
% same
%
% Revision 1.19  2004/03/24 01:25:49  arno
% output only if gui
%
% Revision 1.18  2003/11/05 18:50:22  arno
% fixing normlen problem
%
% Revision 1.17  2003/10/31 19:00:41  arno
% adding more options
%
% Revision 1.16  2003/10/29 16:37:44  arno
% space in command
%
% Revision 1.15  2003/10/29 16:36:49  arno
% typo last
%
% Revision 1.14  2003/10/29 16:35:29  arno
% moving sphere options
%
% Revision 1.13  2003/10/29 16:29:53  arno
% updating command output
%
% Revision 1.12  2003/10/29 16:12:19  arno
% default sphere for dipfit
%
% Revision 1.11  2003/08/04 21:19:14  arno
% command line call bug
%
% Revision 1.10  2003/08/04 18:59:03  arno
% now plot for scanned dipoles
%
% Revision 1.9  2003/06/12 23:51:09  arno
% put normlen by default
%
% Revision 1.8  2003/06/12 23:42:12  arno
% dipfit dipole localization
%
% Revision 1.7  2003/05/06 17:13:32  arno
% debuging passing parameters
%
% Revision 1.6  2003/03/14 22:40:24  arno
% error fif besa head used for dipfit
%
% Revision 1.5  2003/03/14 00:53:04  arno
% debug for dipfit
%
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
    
    geometry = { [2 1] [2 1] [2 1] [2 1] [2.05 0.23 .75] [2.05 0.23 .75] [2.05 0.23 .75] [2.05 0.23 .75] [2.05 0.23 .75] [2.05 0.23 .75] [2 1] };
    uilist = { { 'style' 'text' 'string' 'Components indices ([]=all avaliable)' } ...
               { 'style' 'edit' 'string' '' } ...
               { 'style' 'text' 'string' 'Use dipoles from (scroll, then click to select)' } ...
               { 'style' 'listbox' 'string' varstr(1:end-1) 'value' 1 } ...
               { 'style' 'text' 'string' 'Background image (click to select)' } ...
               { 'style' 'listbox' 'string' fastif(strcmpi(typedip, 'dipfit'), 'Average MRI', 'average MRI|BESA head') } ...
               { 'style' 'text' 'string' 'Plot dipoles within RV (%) range ([min max])' } ...
               { 'style' 'edit' 'string' '' } ...
               { 'style' 'text' 'string' 'Plot sumary mode' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Plot edges' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Plot dipole''s 2-D projections' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Plot projection lines' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Make all dipoles point out' } ...
               { 'style' 'checkbox' 'string' '' } {} ...
               { 'style' 'text' 'string' 'Normalized dipole length' } ...
               { 'style' 'checkbox' 'string' '' 'value' 1 } {} ...
               { 'style' 'text' 'string' 'Additionnal dipplot() options' } ...
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
    options = { options{:} 'image' fastif(result{ind} == 1, 'mri', 'besa') };
    if ~isempty(result{ind+1}), options = { options{:} 'rvrange' eval(  [ '[' result{ind+1} ']' ] ) }; end;
    if result{ind+2} == 1, options = { options{:} 'summary'   'on' }; end;
    if result{ind+3} == 1, options = { options{:} 'drawedges' 'on' }; end;
    if result{ind+4} == 1, options = { options{:} 'projimg'   'on' }; end;
    if result{ind+5} == 1, options = { options{:} 'projlines' 'on' }; end;
    if result{ind+6} == 1, options = { options{:} 'pointout'  'on' }; end; 
    if result{ind+7} == 1, options = { options{:} 'normlen'   'on' }; end;
    if ~isempty( result{ind+8} ), tmpopt = eval( [ '{' result{ind+8} '}' ] ); options = { options{:} tmpopt{:} }; end;
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

    % components to plot
    % ------------------
    if ~isempty(comps)
        if ~isfield(EEG.dipfit.model, 'component')
            for index = comps(:)'
                EEG.dipfit.model(index).component = index;
            end;
        end;
    else
        % find localized dipoles
        comps = [];
        for index2 = 1:length(EEG.dipfit.model)
            if ~isempty(EEG.dipfit.model(index2).posxyz) & EEG.dipfit.model(index2).posxyz(1) ~= 0
                comps = [ comps index2 ];
                EEG.dipfit.model(index2).component = index2;
            end;
        end;        
    end;
    
    % plotting
    % --------
    if ~isempty(findstr(EEG.dipfit.hdmfile, 'BESA'))
        dipplot(EEG.dipfit.model(comps), 'mri', EEG.dipfit.mrifile, 'sph2spm', sph2spm, options{:});
    else
        dipplot(EEG.dipfit.model(comps), 'mri', EEG.dipfit.mrifile, options{:});
    end;
end;
    
if nargin < 3
    com = sprintf('pop_dipplot( %s,%s);', inputname(1), vararg2str({ typedip comps options{:}}));
end;
return;
