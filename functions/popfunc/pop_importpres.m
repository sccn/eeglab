% pop_importpres() - append Presentation event file information into eeglab()
%
% Usage:
%   >> EEGOUT = pop_importpres( EEGIN, filename );
%   >> EEGOUT = pop_importpres( EEGIN, filename, typefield, ...
%                                  latfield, durfield, align, 'key', 'val', ... );
%
% Inputs:
%   EEGIN          - input dataset
%   filename       - file name
%   typefield      - [string] type field name. Default is 'code'.
%   latfield       - [string] latency field name. Default is 'time'.
%   durfield       - [string] duration field name. Default is 'none'.
%   align          - [integer] alignment with preexisting events
%                    see pop_importevent().
%   'key','val'    - This function calls pop_importevent(). These are
%                    optional arguments for this function (for event 
%                    alignment for instance).
% 
% Outputs:
%   EEGOUT         - data structure
%
% Note: 1) If they are pre-existing events in the input dataset,
%          this function will recalculate the latency of the events
%          in the Presentation file, so that they match the one
%          of the pre-existing events.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 15 March 2002
%
% Note: this function is backward compatible (before the input
%       'durfield' was introduced) and can take 'align' as the 
%       5th paramter (instead of 6th parameter).
%
% See also: eeglab(), pop_importevent()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 13 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.15  2004/03/19 00:48:16  arno
% propagate varargin
%
% Revision 1.14  2004/01/29 01:51:16  arno
% debug cancel button
%
% Revision 1.13  2004/01/29 01:39:34  arno
% debug history
%
% Revision 1.12  2004/01/29 01:35:57  arno
% allow to import any field for type and latency
%
% Revision 1.11  2003/12/11 19:52:42  arno
% allowing different alignments
%
% Revision 1.10  2003/11/19 19:28:29  arno
% allow to import different fields for type
%
% Revision 1.9  2003/11/19 18:28:57  arno
% debug code and event type
%
% Revision 1.8  2003/11/07 02:10:39  arno
% debuging decoding of fields
%
% Revision 1.7  2003/11/07 01:45:11  arno
% reading fields automatically
%
% Revision 1.6  2003/11/04 01:11:47  arno
% change default file filter
%
% Revision 1.5  2003/01/24 04:14:46  scott
% header edit -sm
%
% Revision 1.4  2003/01/24 02:14:16  arno
% changing delimiter to 9
%
% Revision 1.3  2002/10/15 17:05:02  arno
% drawnow
%
% Revision 1.2  2002/08/07 17:40:24  arno
% header
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

function [EEG, command] = pop_importpres(EEG, filename, typefield, latfield, durfield, align, varargin); 
command = '';

if nargin < 1 
    help pop_importpres;
    return
end;

% decode input (and backward compatibility)
% -----------------------------------------
if nargin < 5
    durfield = '';
end;
if nargin >= 5 & ~isstr(durfield)
    if nargin >= 6
        varargin = { align varargin{:} };
    end;
    align = durfield;
    durfield = '';
else
    if nargin < 6
        align = 0;
    end;
end;

if nargin < 2 
	% ask user
	[filename, filepath] = uigetfile('*.log;*.LOG', 'Choose a Presentation file -- pop_importpres()'); 
    drawnow;
	if filename == 0 return; end;
	filename = [filepath filename];
end;

fields = loadtxt(filename, 'delim', 9, 'skipline', -2, 'nlines', 1, 'verbose', 'off');

% finding fields
% --------------
if nargin > 1
    if nargin < 3
        typefield = 'code'; % backward compatibility
        latfield  = 'time';
    end;
    indtype  = strmatch(lower(typefield), lower(fields));
    indlat   = strmatch(lower(latfield) , lower(fields));
    if ~isempty(durfield)
         inddur   = strmatch(lower(durfield) , lower(fields));
    else inddur = 0;
    end;
else
    indtype1   = strmatch('event type', lower(fields));
    indtype2   = strmatch('code', lower(fields));
    indtype    = [ indtype1 indtype2 1];
    indlatency = strmatch('time', lower(fields), 'exact');
    indlatency = [ indlatency 1 ];
    uilist = { { 'style' 'text' 'string' [ 'File field containing event types' 10 '' ] } ...
               { 'style' 'list' 'string' strvcat(fields) 'value' indtype(1)  'listboxtop' indtype(1)} ...
               { 'style' 'text' 'string' [ 'File field containing event latencies' 10 '' ] } ...
               { 'style' 'list' 'string' strvcat(fields) 'value' indlatency(1) 'listboxtop' indlatency(1) } ...
               { 'style' 'text' 'string' [ 'File field containing event durations' 10 '' ] } ...
               { 'style' 'list' 'string' strvcat({ 'None' fields{:} }) 'value' 1 'listboxtop' 1 } ...
               { } { 'style' 'text' 'string' 'Note: scroll lists then click to select field' } };
    uigeom = { [2 1] [2 1] [2 1] 1 1 };
    result = inputgui(uigeom, uilist, 'pophelp(''pop_importpres'')', 'Import presentation file - pop_importpres()', ...
                      [], 'normal', [2 2 2 1 1]);
    if isempty(result), return; end;
    
    indtype = result{1};
    indlat  = result{2};
    inddur  = result{3}-1;
    typefield = fields{indtype};
    latfield  = fields{indlat};
    if inddur ~= 0
         durfield  = fields{inddur};
    else durfield  = '';
    end;
end;
if isempty(indtype)
    error(['Could not detect field ''' typefield ''', try importing the file as ASCII (use delimiter=9 (tab))']);
end;
if isempty(indlat)
    error(['Could not detect field ''' latfield ''', try importing the file as ASCII (use delimiter=9 (tab))']);
end;
disp(['Replacing field ''' typefield ''' by ''type'' for EEGLAB compatibility']);
disp(['Replacing field ''' latfield  ''' by ''latency'' for EEGLAB compatibility']);
fields{indtype} = 'type';
fields{indlat}  = 'latency';
if inddur ~= 0
    fields{inddur}  = 'duration';
end;

% regularizing field names
% ------------------------
for index = 1:length(fields)
    indspace = find(fields{index} == ' ');
    fields{index}(indspace) = '_';
    indparen = find(fields{index} == ')');
    if ~isempty(indparen) & indparen == length(fields{index})
        % remove text for parenthesis
        indparen = find(fields{index} == '(');
        if indparen ~= 1
            disp([ 'Renaming ''' fields{index} ''' to ''' fields{index}(1:indparen-1) ''' for Matlab compatibility' ]);
            fields{index} = fields{index}(1:indparen-1);
        else
            fields{index}(indspace) = '_';
        end;
    else
        fields{index}(indspace) = '_';
        indparen = find(fields{index} == '(');
        fields{index}(indspace) = '_';
    end;
end;

% find if uncertainty is duplicated
% ---------------------------------
induncert  = strmatch('uncertainty', lower(fields), 'exact');
if length(induncert) > 1
    fields{induncert(2)}= 'Uncertainty2';
    disp('Renaming second ''Uncertainty'' field');
end;

% import file
% -----------
if isempty(EEG.event), align = NaN; end;

EEG = pop_importevent(EEG, 'append', 'no', 'event', filename, 'timeunit', 1E-4, 'skipline', -3, ...
                           'delim', 9, 'align', align, 'fields', fields, varargin{:});

command = sprintf('EEG = pop_importpres(%s, %s);', inputname(1), vararg2str({ filename typefield latfield durfield align })); 

return;
