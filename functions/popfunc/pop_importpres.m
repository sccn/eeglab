% pop_importpres() - append Presentation event file information into eeglab()
%
% Usage:
%   >> EEGOUT = pop_importpres( EEGIN, filename );
%
% Inputs:
%   EEGIN          - input dataset
%   filename       - file name
%   typefield      - [string] type field name. Default is 'code'.
% 
% Outputs:
%   EEGOUT         - data structure
%
% Note: 1) If they are pre-existing events in the input dataset,
%          this function will recalculate the latency of the events
%          in the Presentation file, so that they match the one
%          of the pre-existing events.
%       2) This function calls pop_importevent()
%
% Author: Arnaud Delorme, CNL / Salk Institute, 15 March 2002
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

function [EEG, command] = pop_importpres(EEG, filename, typefield); 
command = '';

if nargin < 1 
    help pop_importpres;
    return
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
    end;
    indtype  = strmatch(lower(typefield), lower(fields));
else
    indtype  = strmatch('event type', lower(fields));
    indtype2 = strmatch('code', lower(fields));
    if ~isempty(indtype) & ~isempty(indtype2) 
        typefield = questdlg2(strvcat('What field (column) in the .log file do you','want to use for EEGLAB event type?'), ...
                       'pop_importpres()', ...
                       'code','event type','code');
        indtype = strmatch(lower(res), lower(fields));
    else 
        indtype = [indtype indtype2];
    end;        
end;
if isempty(indtype)
    error(['Could not detect field ''' typefield ''', try importing the file as ASCII (use delimiter=9 (tab))']);
end;
disp('Replacing field ''Event Type'' by ''type'' for EEGLAB compatibility');
indlat  = strmatch('time', lower(fields), 'exact');
if isempty(indlat)
    error('Could not detect field ''Time'', try importing the file as ASCII (use delimitor=9 (tab))');
end;
disp('Replacing field ''Time'' by ''latency'' for EEGLAB compatibility');
fields{indtype} = 'type';
fields{indlat}  = 'latency';

% regularizing field names
% ------------------------
for index = 1:length(fields)
    indspace = find(fields{index} == ' ');
    fields{index}(indspace) = '_';
    indparen = find(fields{index} == ')');
    if indparen == length(fields{index})
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
if isempty(EEG.event), align = NaN; 
else                   align = 0; end;
EEG = pop_importevent(EEG, 'append', 'no', 'event', filename, 'timeunit', 1E-4, 'skipline', -3, ...
                           'delim', 9, 'align', align, 'fields', fields);

command = sprintf('EEG = pop_importpres(%s, ''%s'');', inputname(1), filename); 

return;
